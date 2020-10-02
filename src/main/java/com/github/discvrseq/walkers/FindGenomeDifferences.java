package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFRecordCodec;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * This tool is designed to compare two genomes, writing a VCF with positions where the reference differs.  It iterates each base of the source genome, attempts to lift to the target genome, and if liftover is successful it will compare the reference bases.
 * Any positions with differing reference bases will be written out as a VCF file, in the coordinates of the target genome.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar FindGenomeDifferences \
 *     -s sourceGenome.fasta \
 *     -t targetChain.fasta \
 *     -c chainFile.chain \
 *     -O output.vcf.gz
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool can be used to compare two genomes, creating a VCF of positions where the reference alleles differ.",
        oneLineSummary = "Create VCF summarizing positions where reference alleles differ between genomes",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class FindGenomeDifferences extends GATKTool {
    public static final String ORIGINAL_CONTIG = "OriginalContig";
    public static final String ORIGINAL_START = "OriginalStart";

    @Argument(fullName = "sourceGenome", shortName = "s", doc = "Source reference sequence file", optional = false)
    private String sourceGenomeFileName;

    @Argument(fullName = "targetGenome", shortName = "t", doc = "Target reference sequence file", optional = false)
    private String targetGenomeFileName;

    @Argument(fullName = "chainFile", shortName = "c", doc = "Chain file", optional = false)
    private String chainFileName;

    @Argument(doc="File to which the output VCF should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(fullName = "liftoverMinMatch", shortName = "mm", doc = "LiftOver min match", optional = true)
    private double liftoverMinMatch = 0.95;

    private SAMSequenceDictionary sourceDict;
    private SAMSequenceDictionary targetDict;

    private LiftOver liftOver = null;

    private IndexedFastaSequenceFile sourceFasta;
    private IndexedFastaSequenceFile targetFasta;

    private VCFHeader vcfHeader;
    private SortingCollection<VariantContext> sorter;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsReadable(new File(sourceGenomeFileName));
        IOUtil.assertFileIsReadable(new File(targetGenomeFileName));
        IOUtil.assertFileIsReadable(new File(chainFileName));

        IOUtil.assertFileIsWritable(new File(outFile));

        liftOver = new LiftOver(new File(chainFileName));
        liftOver.setLiftOverMinMatch(liftoverMinMatch);
        liftOver.setShouldLogFailedIntervalsBelowThreshold(false);

        sourceDict = SAMSequenceDictionaryExtractor.extractDictionary(new File(sourceGenomeFileName).toPath());
        targetDict = SAMSequenceDictionaryExtractor.extractDictionary(new File(targetGenomeFileName).toPath());

        try {
            sourceFasta = new IndexedFastaSequenceFile(new File(sourceGenomeFileName));
            targetFasta = new IndexedFastaSequenceFile(new File(targetGenomeFileName));
        }
        catch (FileNotFoundException e)
        {
            throw new UserException.BadInput(e.getMessage(), e);
        }

        prepareVcfHeader();
        initializeSorter();
    }

    @Override
    public void traverse() {
        sourceDict.getSequences().forEach(sr -> {
            logger.info("Starting contig: " + sr.getSequenceName());

            final int length = sr.getSequenceLength();
            final byte[] sourceSeq = sourceFasta.getSequence(sr.getSequenceName()).getBases();

            AtomicInteger passedLiftover = new AtomicInteger(0);
            AtomicInteger refDiffers = new AtomicInteger(0);

            IntStream.range(1, length + 1).forEach(sourcePos1 -> {
                Interval lifted = liftOver.liftOver(new Interval(sr.getSequenceName(), sourcePos1, sourcePos1));
                if (lifted != null)
                {
                    passedLiftover.getAndIncrement();

                    byte sourceBase = sourceSeq[sourcePos1 - 1];
                    byte targetBase = getTargetChar(lifted.getContig(), lifted.getStart());
                    if (lifted.isNegativeStrand())
                    {
                        targetBase = BaseUtils.simpleReverseComplement(new byte[]{targetBase})[0];
                    }

                    if (!BaseUtils.basesAreEqual(sourceBase, targetBase))
                    {
                        VariantContextBuilder vcb = new VariantContextBuilder();

                        vcb.chr(lifted.getContig());
                        vcb.start(lifted.getStart());
                        vcb.stop(lifted.getEnd());
                        vcb.alleles(Arrays.asList(Allele.create(targetBase, true), Allele.create(sourceBase)));
                        vcb.attribute(ORIGINAL_CONTIG, sr.getSequenceName());
                        vcb.attribute(ORIGINAL_START, sourcePos1);

                        refDiffers.getAndIncrement();
                        sorter.add(vcb.make());
                    }
                }
            });

            logger.info("total bases: " + length);
            logger.info("total lifted: " + passedLiftover.get());
            logger.info("total lifted with different reference: " + refDiffers.get());
        });

        sorter.doneAdding();
    }

    @Override
    public Object onTraversalSuccess() {
        writeSortedOutput(vcfHeader, sorter);

        return super.onTraversalSuccess();
    }

    private Map<String, byte[]> targetSeqToBytes = new HashMap<>();

    private byte getTargetChar(String targetName, int pos1)
    {
        if (targetSeqToBytes.get(targetName) == null)
        {
            targetSeqToBytes.put(targetName, targetFasta.getSequence(targetName).getBases());
        }

        return targetSeqToBytes.get(targetName)[pos1 - 1];
    }

    private void initializeSorter() {
        File tmpDir = IOUtil.getDefaultTmpDir();
        if (!tmpDir.exists()) {
            tmpDir.mkdirs();
        }

        sorter = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(vcfHeader, true),
                vcfHeader.getVCFRecordComparator(),
                SAMFileWriterImpl.getDefaultMaxRecordsInRam(), tmpDir.toPath());
    }

    private void prepareVcfHeader(){
        vcfHeader = new VCFHeader();
        vcfHeader.setSequenceDictionary(targetDict);

        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(ORIGINAL_CONTIG, 1, VCFHeaderLineType.String, "The name of the contig/chromosome in the source genome."));
        vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(ORIGINAL_START, 1, VCFHeaderLineType.Integer, "The start position of the variant in the source genome."));
    }

    private void writeSortedOutput(final VCFHeader outputHeader, final SortingCollection<VariantContext> sortedOutput) {
        final Log log = Log.getInstance(FindGenomeDifferences.class);
        final ProgressLogger writeProgress = new ProgressLogger(log, 25000, "wrote", "variants");

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setReferenceDictionary(targetDict);
        builder.setOutputFile(new File(outFile));

        try (VariantContextWriter writer = builder.build()) {
            writer.writeHeader(outputHeader);

            for (final VariantContext variantContext : sortedOutput) {
                writer.add(variantContext);
                writeProgress.record(variantContext.getContig(), variantContext.getStart());
            }
        }

        sortedOutput.cleanup();
    }
}
