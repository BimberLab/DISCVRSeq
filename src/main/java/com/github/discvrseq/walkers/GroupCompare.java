package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.github.discvrseq.util.VariableOutputUtils;
import com.opencsv.ICSVWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import javax.annotation.Nullable;
import java.io.IOException;
import java.io.Serial;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * This tool is designed to assist with the generation of candidate variants, based on a given group of samples (e.g., list of samples with a phenotype of interest). It will
 * subset the input VCF to just sites variable in your group(s) of interest, creating a subset VCF with pruned alleles and annotations. Optionally, you can provide a reference
 * VCF (e.g., VCF with genotypes in a broader population). As part of the subsetting, the tool adds several additional INFO field annotations
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar GroupCompare \
 *     -R currentGenome.fasta \
 *     -V myVCF.vcf \
 *     -RV refVCF.vcf \
 *     -O output.vcf.gz \
 *     -G1 'Sample1'
 *     -G1 'Sample2'
 *     -G1 'Sample3'
 *     -G2 'Sample4'
 *     -G2 'Sample5'
 * </pre>
 *
 * <h3>For more advanced usage, you can supply a set of JEXL expressions that will be used to write a subset of variants to a TSV file.</h3>
 * <pre>
 *  java -jar DISCVRseq.jar GroupCompare \
 *     -R currentGenome.fasta \
 *     -V myVCF.vcf \
 *     -RV refVCF.vcf \
 *     -O output.vcf.gz \
 *     -G1 'Sample1'
 *     -G1 'Sample2'
 *     -G1 'Sample3'
 *     -G2 'Sample4'
 *     -G2 'Sample5'
 *     -select 'G1_AF < 0.1 && G2_AF > 0.25'
 *     -select 'IMPACT == "HIGH"'
 *     -select 'CADD_PH > 25'
 *     -F IMPACT
 *     -F CADD_PH
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This is designed to perform simple variant-sifting from a VCF.",
        oneLineSummary = "Produce a list of candidate variants based on a sample list",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class GroupCompare extends ExtendedMultiVariantWalkerGroupedOnStart {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath outFile = null;

    @Argument(doc="File to which a table of high-impact variants should be written. This will only be used if one or more JEXL select expressions are provided.", fullName = "table-output", shortName = "OT", optional = true)
    public GATKPath tableOutput = null;

    @Argument(doc="Group 1", fullName = "group1", shortName = "G1", optional = false)
    public Set<String> group1 = new LinkedHashSet<>(0);

    @Argument(doc="Group 2", fullName = "group2", shortName = "G2", optional = true)
    public Set<String> group2 = new LinkedHashSet<>(0);

    @Argument(fullName= SelectVariants.SELECT_NAME, doc="A filtering expression in terms of either INFO fields or the VariantContext object). If the expression evaluates to true for a variant, it will output in the tabular output.", optional=true)
    private ArrayList<String> selectExpressions = new ArrayList<>();

    @Argument(fullName= "fields-to-output", shortName = "F", doc="The names of fields to include in the output table. These follow the same conventions as GATK VariantsToTable", optional=true)
    private ArrayList<String> additionalFields = new ArrayList<>();

    private VCFHeader outputHeader;

    private List<VCFInfoHeaderLine> getHeaderLines(String prefix) {
        return Arrays.asList(
                new VCFInfoHeaderLine(prefix + "_AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "The name of the contig/chromosome after liftover."),
                new VCFInfoHeaderLine(prefix + "_HOMVAR", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "The start position of the variant on the lifted contig."),
                new VCFInfoHeaderLine(prefix + "_HET", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "The start position of the variant on the lifted contig."),
                new VCFInfoHeaderLine(prefix + "_HOMREF", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "The end position of the variant on the lifted contig.")
        );
    }

    private final String GROUP1 = "G1";
    private final String GROUP2 = "G2";
    private final String REF = "REF";

    final Set<String> samplesToInclude = new LinkedHashSet<>();
    final Set<String> refSamplesToInclude = new LinkedHashSet<>();

    VariantContextWriter writer;
    ICSVWriter csvWriter;

    private List<VariantContextUtils.JexlVCMatchExp> infoJexls = null;

    private List<String> getDefaultOutputFields() {
        List<String> fields = new ArrayList<>(Arrays.asList(
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "IMPACT",
                "G1_AF",
                "G1_HOMREF",
                "G1_HET",
                "G1_HOMVAR"
        ));

        if (group2 != null && !group2.isEmpty()) {
            fields.addAll(Arrays.asList(
                    "G2_AF",
                    "G2_HOMREF",
                    "G2_HET",
                    "G2_HOMVAR"
            ));
        }

        if (getVariantConcordanceScoreArgumentCollection().referenceVariants != null) {
            fields.addAll(Arrays.asList(
                    "REF_AF",
                    "REF_HOMREF",
                    "REF_HET",
                    "REF_HOMVAR"
            ));
        }

        return fields;
    }

    private final List<String> fieldsForTable = new ArrayList<>();

    private final List<String> asFieldsForTable = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        infoJexls = VariantContextUtils.initializeMatchExps(IntStream.range(0, selectExpressions.size()).mapToObj(x -> "Select" + x).toList(), selectExpressions);

        prepareVcfHeader();

        if (tableOutput != null) {
            IOUtil.assertFileIsWritable(tableOutput.toPath().toFile());
            List<String> fields = getDefaultOutputFields();
            fields.addAll(additionalFields);
            fields.forEach(f -> {
                if (VariableOutputUtils.hasGetter(f)) {
                    fieldsForTable.add(f);
                }
                else if (outputHeader.hasInfoLine(f)) {
                    if (outputHeader.getInfoHeaderLine(f).getCountType() == VCFHeaderLineCount.A || outputHeader.getInfoHeaderLine(f).getCountType() == VCFHeaderLineCount.R) {
                        asFieldsForTable.add(f);
                    }
                    else {
                        fieldsForTable.add(f);
                    }
                }
                else {
                    logger.warn("Unknown field: " + f);
                }
            });

            csvWriter = CsvUtils.getTsvWriter(tableOutput.toPath().toFile());
        }
    }

    private void prepareVcfHeader(){
        VCFHeader header = (VCFHeader)getHeaderForFeatures(getVariantConcordanceScoreArgumentCollection().inputVariants);

        Set<VCFHeaderLine> lines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
        lines.addAll(getHeaderLines(GROUP1));
        if (group2 != null && !group2.isEmpty()) {
            lines.addAll(getHeaderLines(GROUP2));
        }
        if (getVariantConcordanceScoreArgumentCollection().referenceVariants != null) {
            lines.addAll(getHeaderLines(REF));
        }

        List<String> sampleNames = new ArrayList<>(group1);
        if (group2 != null) {
            sampleNames.addAll(group2);
        }

        outputHeader = new VCFHeader(lines, sampleNames);
        outputHeader.setSequenceDictionary(getBestAvailableSequenceDictionary());

        samplesToInclude.addAll(group1);
        if (group2 != null) {
            samplesToInclude.addAll(group2);
        }

        if (getVariantConcordanceScoreArgumentCollection().referenceVariants != null) {
            VCFHeader refHeader = (VCFHeader) getHeaderForFeatures(getVariantConcordanceScoreArgumentCollection().referenceVariants);
            refSamplesToInclude.addAll(refHeader.getGenotypeSamples());
            refSamplesToInclude.removeAll(samplesToInclude);
            logger.info("A total of " + refSamplesToInclude.size() + " samples will be selected from the reference VCF");
        }

        List<String> missingSamples = samplesToInclude.stream().filter(Predicate.not(header.getGenotypeSamples()::contains)).toList();
        if (!missingSamples.isEmpty()) {
            throw new GATKException("The following samples were requested but not in the input VCF: " + StringUtils.join(missingSamples, ", "));
        }

        writer = createVCFWriter(outFile);
        writer.writeHeader(outputHeader);
    }

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        Map<FeatureInput<VariantContext>, List<VariantContext>> variants = groupVariantsByFeatureInput(variantContexts);

        VariantContext refVC = null;
        if (variants.get(getVariantConcordanceScoreArgumentCollection().referenceVariants) != null) {
            List<VariantContext> refVCs = variants.get(getVariantConcordanceScoreArgumentCollection().referenceVariants).stream().filter(Predicate.not(VariantContext::isFiltered)).map(x -> x.subContextFromSamples(refSamplesToInclude)).toList();
            refVC = refVCs.size() == 1 ? refVCs.get(0) : GATKVariantContextUtils.simpleMerge(refVCs, Collections.singletonList(getVariantConcordanceScoreArgumentCollection().referenceVariants.getName()), 1, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, true);

            variants.remove(getVariantConcordanceScoreArgumentCollection().referenceVariants);
        }


        for (FeatureInput<VariantContext> fi : variants.keySet()) {
            for (VariantContext vc : variants.get(fi)) {
                VariantContext vc1 = vc.subContextFromSamples(samplesToInclude, true);
                if (!vc1.isVariant() || vc1.isFiltered()) {
                    continue;
                }

                outputVC(vc1, vc, refVC);
            }
        }
    }

    private void outputVC(VariantContext newVC, VariantContext origVC, @Nullable VariantContext refVC) {
        VariantContextBuilder vcb = new VariantContextBuilder(newVC);

        // NOTE: use the original VC to make sure allele order is consistent:
        vcb.attributes(subsetAttributes(newVC, origVC));

        // Add counts for group(s):
        appendCounts(newVC.subContextFromSamples(group1, true), GROUP1, vcb);

        if (group2 != null && !group2.isEmpty()) {
            appendCounts(newVC.subContextFromSamples(group2, true), GROUP2, vcb);
        }

        // Add reference:
        if (refVC != null) {
            appendCounts(refVC, REF, vcb);
        }

        VariantContext toWrite = vcb.make();
        writer.add(toWrite);

        if (csvWriter != null && passesJexlFilters(toWrite)) {
            writeVariantToTable(toWrite);
        }
    }

    private boolean hasWrittenHeader = false;

    private void writeVariantToTable(VariantContext vc) {
        if (!hasWrittenHeader) {
            List<String> fields = new ArrayList<>(fieldsForTable);
            fields.addAll(asFieldsForTable);
            csvWriter.writeNext(fields.toArray(new String[0]));

            hasWrittenHeader = true;
        }

        VariableOutputUtils.extractFields(vc, outputHeader, fieldsForTable, asFieldsForTable, false, true).forEach(line -> {
            csvWriter.writeNext(line.toArray(new String[0]));
        });
    }

    private void appendCounts(VariantContext subsetVC, String prefix, VariantContextBuilder vcb) {
        Map<String, Object> counts = VariantContextUtils.calculateChromosomeCounts(subsetVC, new HashMap<>(), true, Collections.emptySet());
        vcb.attribute(prefix + "_" + VCFConstants.ALLELE_FREQUENCY_KEY, counts.get(VCFConstants.ALLELE_FREQUENCY_KEY));

        vcb.attribute(prefix + "_HOMREF", subsetVC.getGenotypes().stream().filter(Predicate.not(Genotype::isFiltered)).filter(Genotype::isCalled).filter(Genotype::isHomRef).count());
        vcb.attribute(prefix + "_HET", subsetVC.getGenotypes().stream().filter(Predicate.not(Genotype::isFiltered)).filter(Genotype::isCalled).filter(Genotype::isHet).count());
        vcb.attribute(prefix + "_HOMVAR", subsetVC.getGenotypes().stream().filter(Predicate.not(Genotype::isFiltered)).filter(Genotype::isCalled).filter(Genotype::isHomVar).count());
    }

    private Map<String, Object> subsetAttributes(VariantContext newVC, VariantContext origVC) {
        Set<String> keys = origVC.getAttributes().keySet();
        Map<String, Object> attrs = new HashMap<>();
        for (final String key : keys) {
            final VCFInfoHeaderLine headerLine = outputHeader.getInfoHeaderLine(key);
            final int[] relevantIndices = newVC.getAlleles().stream().mapToInt(a -> origVC.getAlleles().indexOf(a)).toArray();
            if (headerLine.getCountType() == VCFHeaderLineCount.A || headerLine.getCountType() == VCFHeaderLineCount.R) {
                attrs.put(key, ReferenceConfidenceVariantContextMerger.generateAnnotationValueVector(headerLine.getCountType(),
                        VariantContextGetters.attributeToList(origVC.getAttribute(key)), relevantIndices));
            }
            else {
                attrs.put(key, origVC.getAttribute(key));
            }
        }

        return attrs;
    }

    @Override
    public Object onTraversalSuccess() {
        writer.close();

        try {
            if (csvWriter != null) {
                csvWriter.close();
            }
        }
        catch (IOException e) {
            logger.error("Unable to close CSV Writer", e);
        }

        return null;
    }

    @Override
    protected GroupCompare.GroupCompareArgumentCollection getMultiVariantInputArgumentCollection() {
        return new GroupCompare.GroupCompareArgumentCollection();
    }

    private GroupCompare.GroupCompareArgumentCollection getVariantConcordanceScoreArgumentCollection() {
        return (GroupCompare.GroupCompareArgumentCollection)multiVariantInputArgumentCollection;
    }

    private static final class GroupCompareArgumentCollection extends MultiVariantInputArgumentCollection {
        @Serial
        private static final long serialVersionUID = 1L;

        @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants")
        public FeatureInput<VariantContext> inputVariants;

        @Argument(doc="Reference Variants", fullName = "referenceVariants", shortName = "RV", optional = true)
        public FeatureInput<VariantContext> referenceVariants = null;

        @Override
        public List<GATKPath> getDrivingVariantPaths() {
            List<GATKPath> ret = new ArrayList<>();
            ret.add(inputVariants);
            if (referenceVariants != null) {
                ret.add(referenceVariants);
            }

            return ret;
        }
    }

    private boolean passesJexlFilters(final VariantContext vc){
        if (infoJexls.isEmpty()){
            return true;
        }

        try {
            for (VariantContextUtils.JexlVCMatchExp jexl : infoJexls) {
                if (VariantContextUtils.match(vc, jexl)){
                    return true;
                }
            }
        } catch (IllegalArgumentException e) {
            throw new UserException(e.getMessage() +
                    "\nSee https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions for documentation on using JEXL in GATK", e);
        }

        return false;
    }
}
