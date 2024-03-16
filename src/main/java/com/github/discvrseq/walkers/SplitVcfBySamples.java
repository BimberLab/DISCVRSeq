package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.google.common.collect.Lists;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.RFC4180Parser;
import com.opencsv.RFC4180ParserBuilder;
import com.opencsv.exceptions.CsvValidationException;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * This takes an input VCF and subsets it into new VCFs where each contains a subset of the original samples. You can either provide samplesPerVcf, in which case it will be subset based on # of samples, or you can provide a
 * TSV file listing the output VCF and sample(s) to add to each.
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

    @Argument(doc="The max number of samples to write per VCF", fullName = "samplesPerVcf", optional = true)
    public Integer samplesPerVcf = null;

    @Argument(doc="If the final VCF in the split has fewer than this number of samples, it will be merged with the second to last VCF", fullName = "minAllowableInFinalVcf", optional = true)
    public Integer minAllowableInFinalVcf = null;

    @Argument(doc="If selected, any site in a subset VCF lacking at least one genotype with a variant will be discarded", fullName = "discardNonVariantSites", optional = true)
    public boolean discardNonVariantSites = false;

    /**
     * When this flag is enabled, all alternate alleles that are not present in the (output) samples will be removed.
     * Note that this even extends to biallelic SNPs - if the alternate allele is not present in any sample, it will be
     * removed and the record will contain a '.' in the ALT column. Note also that sites-only VCFs, by definition, do
     * not include the alternate allele in any genotype calls.
     */
    @Argument(fullName="remove-unused-alternates",
            doc="Remove alternate alleles not present in any genotypes", optional=true)
    public boolean removeUnusedAlternates = false;

    @Argument(doc="This is a TSV file with two columns and no header, where column 1 is the filepath for an output VCF. Column 2 is a sample ID to write to this file. The file can contain multiple rows per output file. The purpose is to supply a list of samples->file, such that this tool can read the input VCF once, and write multiple output VCFs at the same time.", fullName = "sample-mapping-file", optional = true)
    public GATKPath sampleMappingFile = null;

    List<List<String>> batches = new ArrayList<>();
    List<VariantContextWriter> writers = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);
        IOUtil.assertDirectoryIsWritable(outFile.toPath());

        if (samplesPerVcf == null && sampleMappingFile == null) {
            throw new GATKException("Must provide either samplesPerVcf or sampleMappingFile");
        }

        VCFHeader header = getHeaderForVariants();
        List<String> samples = header.getSampleNamesInOrder();

        if (sampleMappingFile != null)
        {
            prepareOutputsForSampleFile(samples, header);
        }
        else
        {
            prepareOutputsForBatches(samples, header);
        }
    }

    private void prepareOutputsForSampleFile(List<String> samples, VCFHeader header) {
        IOUtil.assertFileIsReadable(sampleMappingFile.toPath());

        Map<File, Set<String>> fileToSamples = new HashMap<>();
        RFC4180Parser rfc4180Parser = new RFC4180ParserBuilder().withSeparator('\t').build();
        try (CSVReader reader = new CSVReaderBuilder(IOUtil.openFileForBufferedUtf8Reading(sampleMappingFile.toPath().toFile())).withCSVParser(rfc4180Parser).build()) {
            String[] line;
            while ((line = reader.readNext()) != null) {
                if (line.length != 2) {
                    throw new GATKException("Each sampleMapping file line should only have two elements, found: " + StringUtils.join(line, "/"));
                }

                File f = new File(line[0]);
                IOUtil.assertFileIsWritable(f);

                if (!fileToSamples.containsKey(f)) {
                    fileToSamples.put(f, new TreeSet<>());
                }

                if (!samples.contains(line[1])) {
                    throw new GATKException("Unknown sample: " + line[1]);
                }

                fileToSamples.get(f).add(line[1]);

            }
        }
        catch (IOException | CsvValidationException e) {
            throw new GATKException("Unable to parse sampleMappingFile", e);
        }

        for (File outputFile : fileToSamples.keySet()) {
            VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputFile).setReferenceDictionary(getReferenceDictionary()).build();
            writers.add(writer);
            batches.add(new ArrayList<>(fileToSamples.get(outputFile)));

            VCFHeader outputHeader = new VCFHeader(header.getMetaDataInInputOrder(), fileToSamples.get(outputFile));
            writer.writeHeader(outputHeader);
        }
    }

    private void prepareOutputsForBatches(List<String> samples, VCFHeader header) {
        batches.addAll(Lists.partition(new ArrayList<>(samples), samplesPerVcf));
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
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
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

            VCFHeader outputHeader = new VCFHeader(headerLines, new ArrayList<>(batch));
            writer.writeHeader(outputHeader);
        }
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        int idx = 0;
        for (List<String> batch : batches) {
            final VariantContextWriter writer = writers.get(idx);
            idx++;

            final VariantContext sub = variant.subContextFromSamples(new LinkedHashSet<>(batch), removeUnusedAlternates);
            if (discardNonVariantSites) {
                if (sub.getCalledChrCount() == 0) {
                    continue;
                }
                else if (sub.getGenotypes().stream().noneMatch(g -> !g.isFiltered() && g.isCalled() && !g.isHomRef())) {
                    continue;
                }
            }

            writer.add(sub);
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
