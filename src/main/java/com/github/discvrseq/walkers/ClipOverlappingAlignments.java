package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Feature;
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
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * This tool is designed to soft clip any alignments that start or end in the provided set of intervals. It was originally created to clip reads that overlap with amplification primer binding sites.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar ClipOverlappingAlignments \
 *     -R genome.fasta \
 *     -I input.bam \
 *     --clipIntervals blacklist.bed \
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

    @Argument(doc="File to which variants should be written", fullName = "clipIntervals", shortName = "CI", optional = false)
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
    }

    private long totalOverlappingStart = 0;
    private long totalOverlappingEnd = 0;
    private long totalAlignments = 0;

    @SuppressWarnings("uncheccked")
    private FeatureInput<BEDFeature> addFeatures() {
        return (FeatureInput<BEDFeature>)addFeatureInputsAfterInitialization(bedFile.getPath(), "overlapRegions", BEDFeature.class, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        SAMRecord rec = read.convertToSAMRecord(getHeaderForReads());
        Cigar origCigar = rec.getCigar();

        totalAlignments++;

        if (rec.getReadUnmappedFlag()) {
            writer.addAlignment(rec);
            return;
        }

        featureContext.getValues(bedFeatures).forEach(feat -> {
            //Alignment start within region
            if (read.getSoftStart() >= feat.getStart() && read.getSoftStart() <= feat.getEnd())
            {
                int numBasesToClip = rec.getReadPositionAtReferencePosition(feat.getEnd());
                rec.setCigar(softClipStartOfRead(numBasesToClip, rec.getCigar()));
                totalOverlappingStart++;

                validateCigarChange(rec, origCigar);

                addSummary(feat, numBasesToClip, true);
            }

            //Alignment end within region
            if (read.getSoftEnd() >= feat.getStart() && read.getSoftEnd() <= feat.getEnd())
            {
                int readEnd = rec.getReadPositionAtReferencePosition(feat.getStart()) - 1;
                int numBasesToClip = rec.getAlignmentEnd() - feat.getStart();

                rec.setCigar(new Cigar(CigarUtil.softClipEndOfRead(readEnd, rec.getCigar().getCigarElements())));
                totalOverlappingEnd++;

                validateCigarChange(rec, origCigar);


                addSummary(feat, numBasesToClip, false);
            }
        });

        writer.addAlignment(rec);
    }

    //NOTE: eventually remote this, once it has been run on more diverse CIGAR inputs
    private void validateCigarChange(SAMRecord rec, Cigar origCigar) {
        if (rec.getReadBases().length != rec.getCigar().getReadLength()) {
            throw new GATKException("CIGAR Length Not Equal To Read Length, orig: " + origCigar.toString() + ", new: " + rec.getCigar().toString());
        }

        if (origCigar.getReadLength() != rec.getCigar().getReadLength()) {
            throw new GATKException("CIGAR Length Not Equal, orig: " + origCigar.toString() + ", new: " + rec.getCigar().toString());
        }
    }

    private Cigar softClipStartOfRead(int basesToClip, Cigar oldCigar) {
        List<CigarElement> newCigar = new LinkedList<>();
        newCigar.add(new CigarElement(basesToClip, CigarOperator.S));

        int cigarBasesConsumed = 0;
        for (CigarElement c : oldCigar) {
            final CigarOperator op = c.getOperator();
            final int basesConsumed = op.consumesReadBases() ? c.getLength() : 0;
            final int basesConsumedWithOperator = cigarBasesConsumed + basesConsumed;

            //we are past the soft-clip, just write CIGAR:
            if (cigarBasesConsumed >= basesToClip) {
                newCigar.add(c);
            }
            //we are within the soft clip and this element doesnt span the soft-clip region
            else if (basesConsumedWithOperator <= basesToClip) {
                //Do nothing, allow counting of bases in cigarBasesConsumed
            }
            //we are within the soft clip and need to break apart this element:
            else {
                int adjustedLength = basesConsumed - (basesToClip - cigarBasesConsumed);
                newCigar.add(new CigarElement(adjustedLength, op));
            }

            cigarBasesConsumed += basesConsumedWithOperator;
        }

        return new Cigar(newCigar);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Total alignments written: " + totalAlignments);
        logger.info("Total alignments soft clipped at start: " + totalOverlappingStart);
        logger.info("Total alignments soft clipped at end: " + totalOverlappingEnd);

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
