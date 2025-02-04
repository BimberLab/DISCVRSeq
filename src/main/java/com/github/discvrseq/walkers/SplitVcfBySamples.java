package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.google.common.collect.Lists;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.RFC4180Parser;
import com.opencsv.RFC4180ParserBuilder;
import com.opencsv.exceptions.CsvValidationException;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.list.UnmodifiableList;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

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

    @Argument(fullName="recalculate-ac", doc="This will recalculate the AC, AF, and AN values after subsetting. See also --keep-original-ac", optional=true)
    public boolean recalculateChrCounts = false;

    @Argument(fullName="keep-original-ac", doc="Store the original AC, AF, and AN values after subsetting", optional=true)
    public boolean keepOriginalChrCounts = false;

    @Argument(fullName="original-ac-suffix", doc="If --keep-original-ac is selected, the original AC, AF, and AN values will be stored, but with this suffix (e.g., a suffix of .Orig would result in AF -> AF.Orig)", optional=true)
    public String originalChrCountsSuffix = ".Orig";

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
            batches.add(new UnmodifiableList<>(fileToSamples.get(outputFile).stream().toList()));

            VCFHeader outputHeader = new VCFHeader(getHeaderLines(header), fileToSamples.get(outputFile));
            writer.writeHeader(outputHeader);
        }
    }

    private Set<VCFHeaderLine> getHeaderLines(VCFHeader header) {
        Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(header.getMetaDataInInputOrder());

        if (recalculateChrCounts) {
            if (header.getInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY) == null) {
                headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
            }

            if (header.getInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY) == null) {
                headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
            }

            if (header.getInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY) == null) {
                headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
            }

            if (header.getInfoHeaderLine(VCFConstants.DEPTH_KEY) == null) {
                headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
            }
        }

        if (recalculateChrCounts && keepOriginalChrCounts) {
            headerLines.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY + originalChrCountsSuffix, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Original AC"));
            headerLines.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY + originalChrCountsSuffix, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Original AF"));
            headerLines.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY + originalChrCountsSuffix, 1, VCFHeaderLineType.Integer, "Original AN"));
        }

        return headerLines;
    }

    private void prepareOutputsForBatches(List<String> samples, VCFHeader header) {
        Lists.partition(new ArrayList<>(samples), samplesPerVcf).forEach(l -> {
            batches.add(List.copyOf(l));  // this returns an unmodifiable list
        });

        if (batches.size() == 1) {
            throw new GATKException("This split would result in one output VCF, aborting");
        }

        if (minAllowableInFinalVcf != null) {
            int lastIdx = batches.size() - 1;
            if (batches.get(lastIdx).size() < minAllowableInFinalVcf) {
                List<String> toUpdate = new ArrayList<>(batches.get(lastIdx - 1));
                toUpdate.addAll(batches.get(lastIdx));

                batches.set(lastIdx - 1, new UnmodifiableList<>(toUpdate));
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

            VCFHeader outputHeader = new VCFHeader(getHeaderLines(header), new ArrayList<>(batch));
            writer.writeHeader(outputHeader);
        }
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        int idx = 0;
        for (List<String> sampleNames : batches) {
            final VariantContextWriter writer = writers.get(idx);
            idx++;

            final VariantContext sub = subsetGenotypesBySampleNames(variant, new TreeSet<>(sampleNames), removeUnusedAlternates);
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

    // NOTE: this is basically copied from SelectVariants. This is necessary to ensure we also subset INFO fields correctly
    private VariantContext subsetGenotypesBySampleNames(final VariantContext vc, final SortedSet<String> samples, final boolean removeUnusedAlternates) {
        // If no subsetting of samples or alleles happened, exit now
        if (!removeUnusedAlternates && samples.size() == vc.getNSamples()) {
            return vc;
        }

        // strip out the alternate alleles that aren't being used
        final VariantContext sub = vc.subContextFromSamples(samples, removeUnusedAlternates);

        GenotypesContext newGC;
        if (sub.getNAlleles() != vc.getNAlleles()) {
            // fix the PL and AD values if sub has fewer alleles than original vc
            final GenotypesContext subGenotypesWithOldAlleles = sub.getGenotypes();  //we need sub for the right samples, but PLs still go with old alleles
            newGC = sub.getNAlleles() == vc.getNAlleles() ? subGenotypesWithOldAlleles :
                    AlleleSubsettingUtils.subsetAlleles(subGenotypesWithOldAlleles, 0, vc.getAlleles(),
                            sub.getAlleles(), null, GenotypeAssignmentMethod.DO_NOT_ASSIGN_GENOTYPES);
        } else {
            newGC = sub.getGenotypes();
        }

        // since the VC has been subset (either by sample or allele), we need to strip out the MLE tags
        final VariantContextBuilder builder = new VariantContextBuilder(sub);
        builder.rmAttributes(Arrays.asList(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        builder.genotypes(newGC);
        addAnnotations(builder, vc, sub.getSampleNames());

        if (!vc.getAlleles().equals(sub.getAlleles())) {
            int[] oldToKeyNoRef = sub.getAlternateAlleles().stream().map(a -> vc.getAlleleIndex(a) - 1).mapToInt(Integer::intValue).toArray();
            int[] oldToKeyWithRef = sub.getAlleles().stream().map(vc::getAlleleIndex).mapToInt(Integer::intValue).toArray();
            for (String attr : builder.getAttributes().keySet()) {
                VCFInfoHeaderLine line = getHeaderForVariants().getInfoHeaderLine(attr);
                if (line == null) {
                    continue;
                }

                if (line.getCountType() == VCFHeaderLineCount.A) {
                    builder.attribute(attr, getReorderedAttributes(sub.getAttribute(attr), oldToKeyNoRef));
                }
                else if (line.getCountType() == VCFHeaderLineCount.R) {
                    builder.attribute(attr, getReorderedAttributes(sub.getAttribute(attr), oldToKeyWithRef));
                }
            }
        }

        final VariantContext subset = builder.make();
        return removeUnusedAlternates ? GATKVariantContextUtils.trimAlleles(subset,true,true) : subset;
    }

    private void addAnnotations(final VariantContextBuilder builder, final VariantContext originalVC, final Set<String> selectedSampleNames) {
        if (recalculateChrCounts && keepOriginalChrCounts) {
            final int[] indexOfOriginalAlleleForNewAllele;
            final List<Allele> newAlleles = builder.getAlleles();
            final int numOriginalAlleles = originalVC.getNAlleles();

            // if the alleles already match up, we can just copy the previous list of counts
            if (numOriginalAlleles == newAlleles.size()) {
                indexOfOriginalAlleleForNewAllele = null;
            }
            // otherwise we need to parse them and select out the correct ones
            else {
                indexOfOriginalAlleleForNewAllele = new int[newAlleles.size() - 1];
                Arrays.fill(indexOfOriginalAlleleForNewAllele, -1);

                // note that we don't care about the reference allele at position 0
                for (int newI = 1; newI < newAlleles.size(); newI++) {
                    final Allele newAlt = newAlleles.get(newI);
                    for (int oldI = 0; oldI < numOriginalAlleles - 1; oldI++) {
                        if (newAlt.equals(originalVC.getAlternateAllele(oldI), false)) {
                            indexOfOriginalAlleleForNewAllele[newI - 1] = oldI;
                            break;
                        }
                    }
                }
            }

            if (originalVC.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY + originalChrCountsSuffix,
                        getReorderedAttributes(originalVC.getAttribute(VCFConstants.ALLELE_COUNT_KEY), indexOfOriginalAlleleForNewAllele));
            }
            if (originalVC.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                builder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY + originalChrCountsSuffix,
                        getReorderedAttributes(originalVC.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), indexOfOriginalAlleleForNewAllele));
            }
            if (originalVC.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
                builder.attribute(VCFConstants.ALLELE_NUMBER_KEY + originalChrCountsSuffix, originalVC.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
            }
        }

        if (recalculateChrCounts) {
            VariantContextUtils.calculateChromosomeCounts(builder, false);
        }

        boolean sawDP = false;
        int depth = 0;
        for (final String sample : selectedSampleNames ) {
            final Genotype g = originalVC.getGenotype(sample);
            if (!g.isFiltered()) {
                if (g.hasDP()) {
                    depth += g.getDP();
                    sawDP = true;
                }
            }
        }

        if (sawDP) {
            builder.attribute(VCFConstants.DEPTH_KEY, depth);
        }
    }

    private Object getReorderedAttributes(final Object attribute, final int[] oldToNewIndexOrdering) {
        if (oldToNewIndexOrdering == null || attribute == null) {
            return attribute;
        }

        // break the original attributes into separate tokens; unfortunately, this means being smart about class types
        final Object[] tokens;
        if (attribute.getClass().isArray()) {
            tokens = (Object[]) attribute;
        } else if (List.class.isAssignableFrom(attribute.getClass())) {
            tokens = ((List) attribute).toArray();
        } else {
            tokens = attribute.toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        }

        Utils.validateArg(Arrays.stream(oldToNewIndexOrdering).allMatch(index -> index < tokens.length), () ->
                "the old attribute has an incorrect number of elements: " + attribute);
        return Arrays.stream(oldToNewIndexOrdering).mapToObj(index -> tokens[index]).collect(Collectors.toList());
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
