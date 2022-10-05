package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;

/**
 * This takes an input VCF and subsets it into new VCFs where each contains a subset of the original samples.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar SplitVcfBySample \
 *     -V myVCF.vcf \
 *     -O /outputFolder \
 *     -samplesPerVcf 100 \
 *     -minAllowableInFinalVcf 20
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This takes an input VCF and subsets it into new VCFs where each contains a subset of the original samples.",
        oneLineSummary = "Splits a VCF into multiple VCFs based on samples",
        programGroup = DiscvrSeqProgramGroup.class
)
public class SplitVcfBySamples extends VariantWalker {
    @Argument(doc="The output directory where VCFs will be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath outFile = null;

    @Argument(doc="The max number of samples to write per VCF", fullName = "samplesPerVcf", optional = false)
    public Integer samplesPerVcf = null;

    @Argument(doc="If the final VCF in the split has fewer than this number of samples, it will be merged with the second to last VCF", fullName = "minAllowableInFinalVcf", optional = true)
    public Integer minAllowableInFinalVcf = null;

    /**
     * When this flag is enabled, all alternate alleles that are not present in the (output) samples will be removed.
     * Note that this even extends to biallelic SNPs - if the alternate allele is not present in any sample, it will be
     * removed and the record will contain a '.' in the ALT column. Note also that sites-only VCFs, by definition, do
     * not include the alternate allele in any genotype calls.
     */
    @Argument(fullName="remove-unused-alternates",
            doc="Remove alternate alleles not present in any genotypes", optional=true)
    public boolean removeUnusedAlternates = false;

    List<List<String>> batches;
    List<VariantContextWriter> writers = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);
        IOUtil.assertDirectoryIsWritable(outFile.toPath());

        VCFHeader header = getHeaderForVariants();
        List<String> samples = header.getSampleNamesInOrder();

        batches = Lists.partition(new ArrayList<>(samples), samplesPerVcf);
        if (batches.size() == 1) {
            throw new GATKException("This split would result in one output VCF, aborting");
        }

        if (minAllowableInFinalVcf != null) {
            int lastIdx = batches.size() - 1;
            if (batches.get(lastIdx).size() < minAllowableInFinalVcf) {
                batches.get(lastIdx - 1).addAll(batches.get(lastIdx));
                batches.remove(lastIdx);
            }
        }

        int idx = 0;
        File inputFile = getDrivingVariantsFeatureInput().toPath().toFile();
        for (List<String> batch : batches) {
            idx++;
            String name = inputFile.getName();
            boolean isGz = name.toLowerCase().endsWith(".gz");
            if (isGz) {
                name = FilenameUtils.getBaseName(name);
            }

            name = FilenameUtils.getBaseName(name);

            File output = new File(outFile.toPath().toFile(), name + "." + idx + "of" + batches.size() + ".vcf" + (isGz ? ".gz" : ""));
            VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(output).setReferenceDictionary(getReferenceDictionary()).build();
            writers.add(writer);

            VCFHeader outputHeader = new VCFHeader(header.getMetaDataInInputOrder(), batch);
            writer.writeHeader(outputHeader);
        }
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        int idx = 0;
        for (List<String> batch : batches) {
            final VariantContextWriter writer = writers.get(idx);
            final VariantContext sub = variant.subContextFromSamples(new LinkedHashSet<>(batch), removeUnusedAlternates);

            writer.add(sub);
            idx++;
        }
    }

    @Override
    public void closeTool() {
        super.closeTool();

        for (VariantContextWriter writer : writers) {
            try {
                writer.close();
            }
            catch (Exception e) {
                logger.error("Error closing VCF writer", e);
            }
        }
    }
}
