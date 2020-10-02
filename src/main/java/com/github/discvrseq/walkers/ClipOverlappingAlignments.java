package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * This tool is designed to soft clip any alignments that start or end in the provided set of intervals. It was originally created to clip reads that overlap with amplification primer binding sites.
 *
 * Note: this tool has not been tested on every possible BAM/CIGAR configuration and it's possible that in certain cases
 * other information stored in the BAM isnt properly updated.  For example, if a given alignment no longer maps to the reference
 * after clipping it is marked unmapped, but the mate is not changed (which can store information about the mate's position). We therefore
 * recommend running Picard's FixMateInformation afterwards.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar ClipOverlappingAlignments \
 *     -R genome.fasta \
 *     -I input.bam \
 *     --clipIntervals primerSites.bed \
 *     --reportFile summary.txt \
 *     -O output.bam
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool can be used to soft-clip alignments overlapping the specified intervals, such as a BED file of primer binding sites",
        oneLineSummary = "Clip alignments overlapping the specified intervals",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class ClipOverlappingAlignments extends ReadWalker {
    @Argument(doc="File to which alignment should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File outFile = null;

    @Argument(doc="BED file specifying the clipping intervals. Any alignments starting or ending within these intervals will be soft-clipped to these borders", fullName = "clipIntervals", shortName = "CI", optional = false)
    public File bedFile = null;

    @Argument(doc="File to which a summary of clipping will be written", fullName = "reportFile", shortName = "rf", optional = true)
    public File reportFile = null;

    private FeatureInput<BEDFeature> bedFeatures;

    private SAMFileWriter writer;

    private Map<String, HitTracker> summary = new LinkedHashMap<>();

    private static class HitTracker {
        final String contig;
        final int start;
        final int end;

        int totalClipped5 = 0;
        int totalClipped3 = 0;

        int totalBasesClipped5 = 0;
        int totalBasesClipped3 = 0;

        public HitTracker(String contig, int start, int end) {
            this.contig = contig;
            this.start = start;
            this.end = end;
        }
    }

    private void addSummary(BEDFeature feat, int totalClipped, boolean isFivePrime) {
        String key = feat.getContig() + "<>" + feat.getStart() + "<>" + feat.getEnd();
        HitTracker tracker = summary.getOrDefault(key, new HitTracker(feat.getContig(), feat.getStart(), feat.getEnd()));
        if (isFivePrime) {
            tracker.totalClipped5++;
            tracker.totalBasesClipped5 += totalClipped;
        }
        else {
            tracker.totalClipped3++;
            tracker.totalBasesClipped3 += totalClipped;
        }

        summary.put(key, tracker);
    }

    @Override
    public void onTraversalStart() {
        IOUtil.assertFileIsWritable(outFile);
        IOUtil.assertFileIsReadable(bedFile);

        writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(getHeaderForReads(), false, outFile);

        File expected = new File(bedFile.getParentFile(), bedFile.getName() + FileExtensions.TRIBBLE_INDEX);
        if (!expected.exists())
        {
            logger.info("Writing index for file: " + bedFile.getPath());
            Index index = IndexFactory.createDynamicIndex(bedFile, new BEDCodec());
            try {
                index.writeBasedOnFeatureFile(bedFile);
            }
            catch (IOException e){
                throw new RuntimeException(e);
            }
        }

        if (features == null) {
            features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, getGenomicsDBOptions());
        }

        bedFeatures = addFeatures();

        features.getFeatureIterator(bedFeatures).forEachRemaining(feat -> {
            if (feat.getStart() >= feat.getEnd()) {
                throw new UserException.BadInput("Improper BED feature, start greater than end: " + feat.getContig() + ":" + feat.getStart() + "-" + feat.getEnd());
            }
        });
    }

    private long totalOverlappingStart = 0;
    private long totalOverlappingEnd = 0;
    private long totalAlignments = 0;
    private long readsDropped = 0;

    @SuppressWarnings("unchecked")
    private FeatureInput<BEDFeature> addFeatures() {
        return (FeatureInput<BEDFeature>)addFeatureInputsAfterInitialization(bedFile.getPath(), "overlapRegions", BEDFeature.class, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        SAMRecord rec = read.convertToSAMRecord(getHeaderForReads());
        Cigar origCigar = rec.getCigar();

        totalAlignments++;
        AtomicBoolean shouldWrite = new AtomicBoolean(true);

        if (rec.getReadUnmappedFlag()) {
            writer.addAlignment(rec);
            return;
        }

        featureContext.getValues(bedFeatures).forEach(feat -> {
            //Alignment start within region
            if (read.getStart() >= feat.getStart() && read.getStart() <= feat.getEnd())
            {
                int newAlignStart = feat.getEnd() + 1;

                //Increment start until it doesnt land in a deletion:
                while (newAlignStart <= read.getEnd() && rec.getReadPositionAtReferencePosition(newAlignStart) == 0) {
                    newAlignStart++;
                }

                if (newAlignStart >= read.getEnd()) {
                    readsDropped++;
                    setUnaligned(rec);
                    if (rec.isSecondaryOrSupplementary()) {
                        shouldWrite.getAndSet(false);
                    }
                    return;
                }

                int numBasesToClip = rec.getReadPositionAtReferencePosition(newAlignStart) - 1;
                rec.setCigar(softClipStartOfRead(numBasesToClip, rec.getCigar()));
                rec.setAttribute("Xc", origCigar.toString());
                if (rec.getCigar().getReferenceLength() == 0) {
                    readsDropped++;
                    setUnaligned(rec);
                    if (rec.isSecondaryOrSupplementary()) {
                        shouldWrite.getAndSet(false);
                    }
                    return;
                }

                rec.setAlignmentStart(newAlignStart);
                totalOverlappingStart++;

                validateCigarChange(rec, origCigar);
                addSummary(feat, numBasesToClip, true);


            }

            //Alignment end within region
            if (read.getEnd() >= feat.getStart() && read.getEnd() <= feat.getEnd())
            {
                int newAlignEnd = feat.getStart() - 1;

                //Decrement until it doesnt land in a deletion:
                while (newAlignEnd >= read.getStart() && rec.getReadPositionAtReferencePosition(newAlignEnd) == 0) {
                    newAlignEnd--;
                }

                if (newAlignEnd <= read.getStart()) {
                    readsDropped++;
                    if (rec.isSecondaryOrSupplementary()) {
                        shouldWrite.getAndSet(false);
                    }
                    setUnaligned(rec);
                    return;
                }

                int readEnd = rec.getReadPositionAtReferencePosition(newAlignEnd);
                int numBasesToClip = rec.getAlignmentEnd() - newAlignEnd;

                try {
                    rec.setCigar(new Cigar(CigarUtil.softClipEndOfRead(readEnd, rec.getCigar().getCigarElements())));
                }
                catch (SAMException e) {
                    logger.error("Problem clipping end of read: " + read.getName() + ".  CIGAR: " + rec.getCigar().toString() + ", " + readEnd);
                    throw e;
                }

                rec.setAttribute("Xc", origCigar.toString());
                if (rec.getCigar().getReferenceLength() == 0) {
                    readsDropped++;
                    if (rec.isSecondaryOrSupplementary()) {
                        shouldWrite.getAndSet(false);
                    }
                    setUnaligned(rec);
                    return;
                }
                totalOverlappingEnd++;

                validateCigarChange(rec, origCigar);

                addSummary(feat, numBasesToClip, false);
            }
        });

        if (shouldWrite.get()) {
            List<SAMValidationError> errors = rec.isValid();
            if (errors != null && !errors.isEmpty()) {
                for (SAMValidationError e : errors) {
                    logger.error(e.getMessage());
                }

                throw new GATKException("Invalid SAM Record: " + origCigar.toString() + ", new: " + rec.getCigar().toString() + ", align start: " + rec.getAlignmentStart() + ", name: " + rec.getReadName());
            }
            writer.addAlignment(rec);
        }
    }

    private SAMRecord setUnaligned(SAMRecord rec) {
        rec.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        rec.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        rec.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        rec.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
        rec.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        rec.setReadUnmappedFlag(true);

        return rec;
    }

    //NOTE: eventually remote this, once it has been run on more diverse CIGAR inputs
    private void validateCigarChange(SAMRecord rec, Cigar origCigar) {
        if (rec.getReadBases().length != rec.getCigar().getReadLength()) {
            throw new GATKException("CIGAR Read Length Not Equal To Actual Read Length, orig: " + origCigar.toString() + ", new: " + rec.getCigar().toString() + ", align start: " + rec.getAlignmentStart() + ", name: " + rec.getReadName());
        }

        if (origCigar.getReadLength() != rec.getCigar().getReadLength()) {
            throw new GATKException("CIGAR Length Not Equal, orig: " + origCigar.toString() + ", new: " + rec.getCigar().toString() +", align start: " + rec.getAlignmentStart() + ", name: " + rec.getReadName());
        }

        int s = rec.getReadPositionAtReferencePosition(rec.getAlignmentStart());
        if (s == 0) {
            throw new GATKException("Read alignment start doesnt correspond to read bases: " + origCigar.toString() + ", new: " + rec.getCigar().toString() + ", align start: " + rec.getAlignmentStart() + ", name: " + rec.getReadName());
        }

        int e = rec.getReadPositionAtReferencePosition(rec.getAlignmentEnd());
        if (e == 0) {
            throw new GATKException("Read alignment end doesnt correspond to read bases: " + origCigar.toString() + ", new: " + rec.getCigar().toString() + ", align start: " + rec.getAlignmentStart() + ", name: " + rec.getReadName());
        }

        if (rec.getCigar().getReferenceLength() == 0) {
            throw new GATKException("Read has zero reference length: " + origCigar.toString() + ", new: " + rec.getCigar().toString() + ", align start: " + rec.getAlignmentStart() + ", name: " + rec.getReadName());
        }
    }

    private Cigar softClipStartOfRead(int basesToClip, Cigar oldCigar) {
        List<CigarElement> newCigar = new LinkedList<>();
        newCigar.add(new CigarElement(basesToClip, CigarOperator.S));

        int totalReadBasesConsumed = 0;
        boolean outsidePadding = false;
        for (CigarElement c : oldCigar) {
            final CigarOperator op = c.getOperator();
            final int basesConsumedByOperator = op.consumesReadBases() ? c.getLength() : 0;
            final int basesConsumedWithCurrentOperator = totalReadBasesConsumed + basesConsumedByOperator;

            //First non-padded element cannot be a deletion
            if (!outsidePadding && op == CigarOperator.D) {
                continue;
            }

            //we are past the soft-clip, just write CIGAR:
            if (totalReadBasesConsumed >= basesToClip) {
                addOrMergeElement(newCigar, c);
                if (!op.isPadding()) {
                    outsidePadding = true;
                }
            }
            //we are within the soft clip and this element doesnt span the soft-clip region
            else if (basesConsumedWithCurrentOperator <= basesToClip) {
                //Do nothing, allow counting of bases in totalReadBasesConsumed
            }
            //we are within the soft clip and need to break apart this element:
            else {
                int adjustedLength = basesConsumedByOperator - (basesToClip - totalReadBasesConsumed);
                addOrMergeElement(newCigar, new CigarElement(adjustedLength, op));
                if (!op.isPadding()) {
                    outsidePadding = true;
                }
            }

            totalReadBasesConsumed = basesConsumedWithCurrentOperator;
        }

        return new Cigar(newCigar);
    }

    private void addOrMergeElement(List<CigarElement> list, CigarElement toAdd) {
        if (list.isEmpty()) {
            list.add(toAdd);
            return;
        }

        int lastIdx = list.size() - 1;
        CigarElement lastEl = list.get(lastIdx);
        if (lastEl.getOperator() == toAdd.getOperator()) {
            list.set(lastIdx, new CigarElement(toAdd.getLength() + lastEl.getLength(), toAdd.getOperator()));
        }
        else {
            list.add(toAdd);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Total alignments written: " + totalAlignments);
        logger.info("Total alignments soft clipped at start: " + totalOverlappingStart);
        logger.info("Total alignments soft clipped at end: " + totalOverlappingEnd);
        logger.info("Total alignments dropped due to lack of reference coverage after clipping: " + readsDropped);

        if (reportFile != null) {
            try (CSVWriter csvWriter = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(reportFile), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
                csvWriter.writeNext(new String[]{"FeatureContig", "FeatureStart", "FeatureEnd","TotalAlignmentsClippedAt5Prime", "Total5PrimeBasesClipped","TotalAlignmentsClippedAt3Prime", "Total3PrimeBasesClipped"});
                for (String key : summary.keySet()) {
                    HitTracker tracker = summary.get(key);
                    csvWriter.writeNext(new String[]{tracker.contig, String.valueOf(tracker.start), String.valueOf(tracker.end), String.valueOf(tracker.totalClipped5), String.valueOf(tracker.totalBasesClipped5), String.valueOf(tracker.totalClipped3), String.valueOf(tracker.totalBasesClipped3)});
                }
            }
             catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
             }
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        if (writer != null)
            writer.close();
    }
}
