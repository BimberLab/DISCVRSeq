package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import org.apache.logging.log4j.core.util.Assert;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

/**
 *
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar ClipOverlappingAlignments \
 *     -R currentGenome.fasta \
 *     -I input.bam \
 *     --clipIntervals blacklist.bed \
 *     -O output.bam
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class ClipOverlappingAlignments extends ReadWalker {
    @Argument(doc="File to which alignment should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File outFile = null;

    @Argument(doc="File to which variants should be written", fullName = "clipIntervals", shortName = "CI", optional = false)
    public File bedFile = null;

    private FeatureInput<BEDFeature> bedFeatures;

    private SAMFileWriter writer;

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

        bedFeatures = (FeatureInput<BEDFeature>)addFeatureInputsAfterInitialization(bedFile.getPath(), "overlapRegions", BEDFeature.class, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    private long totalOverlappingStart = 0;
    private long totalOverlappingEnd = 0;
    private long totalAlignments = 0;

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        SAMRecord rec = read.convertToSAMRecord(getHeaderForReads());
        Cigar orig = rec.getCigar();

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
            }

            //Alignment end within region
            if (read.getSoftEnd() >= feat.getStart() && read.getSoftEnd() <= feat.getEnd())
            {
                int readEnd = rec.getReadPositionAtReferencePosition(feat.getStart()) - 1;
                rec.setCigar(new Cigar(CigarUtil.softClipEndOfRead(readEnd, rec.getCigar().getCigarElements())));
                totalOverlappingEnd++;
            }
        });

        if (rec.getReadString().length() != rec.getCigar().getReadLength()) {
            throw new GATKException("CIGAR Length Not Equal To Read Length, orig: " + orig.toString() + ", new: " + rec.getCigar().toString());
        }

        if (orig.getReadLength() != rec.getCigar().getReadLength()) {
            throw new GATKException("CIGAR Length Not Equal, orig: " + orig.toString() + ", new: " + rec.getCigar().toString());
        }

        writer.addAlignment(rec);
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

        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        if (writer != null)
            writer.close();
    }
}
