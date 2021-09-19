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
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Merge multiple VCFs to produce one with all unique sites.  All genotypes will be dropped.",
        oneLineSummary = "Merge multiple VCFs to produce one with all unique sites",
        programGroup = VariantManipulationProgramGroup.class
)
public class MergeVariantSites extends MultiVariantWalker {

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
    private boolean EXCLUDE_FILTERED = false;

    private VariantContextWriter vcfWriter;

    @Override
    public void onTraversalStart() {
        //Clone header, without samples
        VCFHeader vcfHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), Collections.emptySet());

        vcfWriter = createVCFWriter(out);
        vcfWriter.writeHeader(vcfHeader);
    }

    private final List<VariantContext> vcsForSite = new ArrayList<>();

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (EXCLUDE_FILTERED && vc.isFiltered()) {
            return;
        }

        if (EXCLUDE_NON_VARIANTS  && !(vc.isPolymorphicInSamples() && !GATKVariantContextUtils.isSpanningDeletionOnly(vc)))
        {
            return;
        }

        if (vcsForSite.isEmpty())
        {
            vcsForSite.add(vc);
            return;
        }

        if (!vc.getContig().equals(vcsForSite.get(0).getContig()) || vc.getStart() != vcsForSite.get(0).getStart()) {
            processAndWriteVariant();
        }

        vcsForSite.add(vc);
    }

    private void processAndWriteVariant()
    {
        VariantContext vc = GATKVariantContextUtils.simpleMerge(vcsForSite, null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_UNCONDITIONAL, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE, true);
        VariantContextBuilder vcb = new VariantContextBuilder(vc.getSource(), vc.getContig(), vc.getStart(), vc.getEnd(), vc.getAlleles());

        vc = vcb.make();
        vcfWriter.add(vc);
        vcsForSite.clear();
    }

    @Override
    public Object onTraversalSuccess() {
        if (!vcsForSite.isEmpty()) {
            processAndWriteVariant();
        }

        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();

        if (vcfWriter != null)
            vcfWriter.close();
    }
}