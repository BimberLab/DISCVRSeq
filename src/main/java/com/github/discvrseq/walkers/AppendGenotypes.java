package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Merge multiple VCFs, potentially with distinct samples.",
        oneLineSummary = "Merge multiple VCFs to produce one with all sites",
        programGroup = VariantManipulationProgramGroup.class
)
public class AppendGenotypes extends VariantWalker {

    @Argument(doc="The files providing genotypes.  The set of samples must be unique across them", fullName = "genotypeVcfs", shortName = "g", optional = false)
    public List<FeatureInput<VariantContext>> genotypeVcfs = new ArrayList<>();

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    private Set<String> samples;

    private VariantContextWriter vcfWriter = null;

    @Override
    public void onTraversalStart() {
        IOUtil.assertFileIsWritable(out);

        List<VCFHeader> headers = features.getAllVariantHeaders();
        Set<String> samples = new HashSet<>();
        for (VCFHeader header : headers) {
            if (CollectionUtils.containsAny(samples, header.getSampleNamesInOrder())) {
                throw new GATKException("All VCFs must have distinct samples.  Overlap in: " + StringUtils.join(CollectionUtils.intersection(samples, header.getSampleNamesInOrder()), ", "));
            }

            samples.addAll(header.getSampleNamesInOrder());
        }

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(headers, true);
        GATKVariantContextUtils.addChromosomeCountsToHeader(headerLines);

        VCFHeader vcfHeader = new VCFHeader(headerLines, samples);

        vcfWriter = createVCFWriter(out);
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        GenotypesContext genotypes = GenotypesContext.create();
        List<VariantContext> toMerge = new ArrayList<>();

        for (FeatureInput<VariantContext> g : genotypeVcfs) {
            List<VariantContext> vcs = featureContext.getValues(g);
            if (vcs.size() > 1) {
                throw new GATKException("More than one variant found for position: " + variant.getContig() + " " + variant.getStart() + " in " + g.getName());
            }
            else if (vcs.isEmpty())
            {
                continue;
            }

            //Check identical Alleles
            VariantContext toAdd = vcs.get(0);
            if (!variant.getAlleles().equals(toAdd.getAlleles())) {
                toMerge.add(toAdd);
            }
            else {
                genotypes.addAll(toAdd.getGenotypes());
            }
        }

        if (!toMerge.isEmpty()) {
            toMerge.add(0, variant);
            variant = GATKVariantContextUtils.simpleMerge(toMerge, null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_UNCONDITIONAL, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE, true);
            genotypes.addAll(variant.getGenotypes());
        }

        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        vcb.genotypes(genotypes);

        //Re-calculate these
        VariantContextUtils.calculateChromosomeCounts(vcb, false);

        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        super.closeTool();

        if (vcfWriter != null)
            vcfWriter.close();
    }
}