package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.milaboratory.core.PairedEndReadsLayout;
import com.milaboratory.core.merger.MergerParameters;
import com.milaboratory.core.merger.MismatchOnlyPairedReadMerger;
import com.milaboratory.core.merger.PairedReadMergingResult;
import com.milaboratory.core.merger.QualityMergingAlgorithm;
import com.milaboratory.core.sequence.NSequenceWithQuality;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;

import java.io.File;

/**
 * This tool accepts a pair of FASTQ files and attempts to merge each read pair into a single read, based on the values for minimumOverlap and minimalIdentity.
 *
 * <h3>Usage example, table output only:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar MergeFastqReads \
 *     -fq1 reads1.fastq.gz \
 *     -fq2 reads2.fastq.gz \
 *     -minLength 100 \
 *     -O merged.fastq.gz
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool iterates a pair of FASTQs and attempts to merge the forward/reverse reads into a single read.",
        oneLineSummary = "Merges paired FASTQ reads into a single read",
        programGroup = DiscvrSeqProgramGroup.class
)
public class MergeFastqReads extends GATKTool {
    @Argument(fullName="fastq1", shortName = "fq1", doc="Input fastq file (optionally gzipped) of forward reads.")
    public File FASTQ;

    @Argument(fullName="fastq2", shortName = "fq2", doc="Input fastq file (optionally gzipped) of reverse reads.")
    public File FASTQ2;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The output FASTQ file")
    private File outputFile1 = null;

    //@Argument(fullName = "qualityMergingAlgorithm", doc = "The algorithm for merging, one of: ")
    private QualityMergingAlgorithm qualityMergingAlgorithm = QualityMergingAlgorithm.MaxSubtraction;

    @Argument(fullName = "minimalOverlap", doc = "The minimum overlap in base pairs")
    private int minimalOverlap = 17;

    @Argument(fullName = "minimalIdentity", doc = "The minimum identity in the overlap region")
    private double minimalIdentity = 0.9;

    @Argument(fullName = "minLength", doc = "The minimum length of the merged reads")
    private double minLength = 0;

    //@Argument(fullName = "identityType", doc = "")
    private String identityType = "Unweighted";


    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsReadable(FASTQ2);
        IOUtil.assertFileIsWritable(outputFile1);
    }

    @Override
    public void traverse() {
        FastqWriterFactory fact = new FastqWriterFactory();
        fact.setUseAsyncIo(true);

        MismatchOnlyPairedReadMerger merger = new MismatchOnlyPairedReadMerger(minimalOverlap, minimalIdentity, MergerParameters.DEFAULT_MAX_QUALITY_VALUE, qualityMergingAlgorithm, PairedEndReadsLayout.Opposite);

        long written = 0L;
        long failed = 0L;
        try (FastqReader reader1 = fileToFastqReader(FASTQ); FastqReader reader2 = FASTQ2 == null ? null : fileToFastqReader(FASTQ2); FastqWriter writer1 = fact.newWriter(outputFile1)) {
            while(reader1.hasNext())
            {
                FastqRecord fq1 = reader1.next();
                FastqRecord fq2 = reader2 == null ? null : reader2.next();

                NSequenceWithQuality ns1 = new NSequenceWithQuality(fq1.getReadString(), fq1.getBaseQualityString());
                NSequenceWithQuality ns2 = new NSequenceWithQuality(fq2.getReadString(), fq2.getBaseQualityString());

                PairedReadMergingResult result = merger.merge(ns1, ns2);
                if (result.isSuccessful()) {
                    NSequenceWithQuality merged = result.getOverlappedSequence();
                    if (merged.getSequence().size() < minLength) {
                        failed++;
                    } else {
                        writer1.write(new FastqRecord(fq1.getReadName(), merged.getSequence().toString(), fq1.getBaseQualityHeader(), merged.getQuality().toString()));
                        written++;
                    }
                }
                else {
                    failed++;
                }
            }
        }

        logger.info("total reads merged: " + written);
        logger.info("failed merge: " + failed);
    }

    private FastqReader fileToFastqReader(final File file) {
        return new FastqReader(file, true);
    }

}
