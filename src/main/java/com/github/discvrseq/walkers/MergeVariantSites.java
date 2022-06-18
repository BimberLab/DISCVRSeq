package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Merge multiple VCFs to produce one with all unique sites.  All genotypes will be dropped.",
        oneLineSummary = "Merge multiple VCFs to produce one with all unique sites",
        programGroup = VariantManipulationProgramGroup.class
)
public class MergeVariantSites extends ExtendedMultiVariantWalkerGroupedOnStart {

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    /**
     * Exclude sites that do not contain any called ALT alleles in the merged callset. The evaluation is made after the
     * merging procedure is complete.
     */
    @Argument(fullName="excludeNonVariants", shortName="env", doc="Exclude sites where no variation is present in the input", optional=true)
    public boolean EXCLUDE_NON_VARIANTS = false;

    /**
     * If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS`
     * in the FILTER field) will be excluded from the output.
     */
    @Argument(fullName="exclude-filtered", doc="Don't include filtered sites", optional=true)
    public boolean EXCLUDE_FILTERED = false;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        //Clone header, without samples
        VCFHeader vcfHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), Collections.emptySet());

        vcfWriter = createVCFWriter(out);
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        List<VariantContext> toMerge = new ArrayList<>(variantContexts);
        if (EXCLUDE_FILTERED) {
            toMerge = toMerge.stream().filter(vc -> !vc.isFiltered()).collect(Collectors.toList());
        }

        if (EXCLUDE_NON_VARIANTS)
        {
            toMerge = toMerge.stream().filter(vc ->  !(vc.isPolymorphicInSamples() && !GATKVariantContextUtils.isSpanningDeletionOnly(vc))).collect(Collectors.toList());
        }

        if (toMerge.isEmpty())
        {
            return;
        }

        List<VariantContext> toMergeNoGenotype = toMerge.stream().map(vc -> new VariantContextBuilder(vc.getSource(), vc.getContig(), vc.getStart(), vc.getEnd(), vc.getAlleles()).make()).collect(Collectors.toList());
        VariantContext vc = GATKVariantContextUtils.simpleMerge(toMergeNoGenotype, null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_UNCONDITIONAL, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE, true);
        vcfWriter.add(vc);
    }

    @Override
    public void closeTool() {
        super.closeTool();

        if (vcfWriter != null)
            vcfWriter.close();
    }
}