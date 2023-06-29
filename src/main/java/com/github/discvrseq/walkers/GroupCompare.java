package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
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
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import javax.annotation.Nullable;
import java.util.*;
import java.util.function.Predicate;

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

    @Argument(doc="Reference Variants", fullName = "referenceVariants", shortName = "RV", optional = true)
    public FeatureInput<VariantContext> referenceVariants = null;

    @Argument(doc="Group 1", fullName = "group1", shortName = "G1", optional = false)
    public Set<String> group1 = new LinkedHashSet<>(0);

    @Argument(doc="Group 2", fullName = "group2", shortName = "G2", optional = true)
    public Set<String> group2 = new LinkedHashSet<>(0);

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

    VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        prepareVcfHeader();
    }

    private void prepareVcfHeader(){
        VCFHeader header = getHeaderForVariants();

        Set<VCFHeaderLine> lines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
        lines.addAll(getHeaderLines(GROUP1));
        if (group2 != null && !group2.isEmpty()) {
            lines.addAll(getHeaderLines(GROUP2));
        }
        if (referenceVariants != null) {
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

        if (referenceVariants != null) {
            VCFHeader refHeader = (VCFHeader) getHeaderForFeatures(referenceVariants);
            List<String> rs = new ArrayList<>(refHeader.getGenotypeSamples());
            rs.removeAll(samplesToInclude);
            logger.info("A total of " + rs.size() + " samples will be selected from the reference VCF");
        }

        List<String> missingSamples = samplesToInclude.stream().filter(Predicate.not(getHeaderForVariants().getGenotypeSamples()::contains)).toList();
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
        if (variants.get(referenceVariants) != null) {
            List<VariantContext> refVCs = variants.get(referenceVariants).stream().filter(Predicate.not(VariantContext::isFiltered)).map(x -> x.subContextFromSamples(samplesToInclude)).toList();
            refVC = refVCs.size() == 1 ? refVCs.get(0) : GATKVariantContextUtils.simpleMerge(refVCs, Collections.singletonList(referenceVariants.getName()), 1, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, true);

            variants.remove(referenceVariants);
        }


        for (FeatureInput<VariantContext> fi : getDrivingVariantsFeatureInputs()) {
            for (VariantContext vc : variants.get(fi)) {
                VariantContext vc1 = vc.subContextFromSamples(samplesToInclude, true);
                if (!vc1.isVariant()) {
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

        writer.add(vcb.make());
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

        return null;
    }
}
