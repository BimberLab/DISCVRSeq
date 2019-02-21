package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqDevProgramGroup;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool accepts a BED file with a list of sites, where sample name is the feature name.  If the VCF has a variant at any of these positions, the sample(s) will have their genotypes converted to no-call.  The BED file can have multiple lines per site if more than one sample needs to be filtered.",
        oneLineSummary = "Filters genotypes based on an input list of samples/sites",
        programGroup = DiscvrSeqDevProgramGroup.class
)
public class GenotypeFilterBySample extends VariantWalker {

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    @Argument(doc = "Genotype Blastlist BED", fullName = "genotypeBlacklist", shortName = "bl", optional = false)
    public FeatureInput<BEDFeature> genotypeBlacklist = null;

    private VariantContextWriter writer = null;

    private long totalFiltered = 0L;

    @Override
    public void onTraversalStart() {
        VCFHeader header = new VCFHeader(getHeaderForVariants());

        writer = createVCFWriter(out);
        writer.writeHeader(header);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        GenotypesContext gc = GenotypesContext.copy(variant.getGenotypes());

        //look for overlapping features in the input
        List<BEDFeature> overlapping = featureContext.getValues(genotypeBlacklist);
        if (!overlapping.isEmpty())
        {
            //parse to find ID(s)
            List<String> sampleIDs = new ArrayList<>();
            for (BEDFeature f : overlapping)
            {
                sampleIDs.add(f.getName());
            }

            for (String sampleId : sampleIDs)
            {
                Genotype g = variant.getGenotype(sampleId);
                if (g != null)
                {
                    GenotypeBuilder gb = new GenotypeBuilder(g);
                    gb.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
                    gc.replace(gb.make());
                    totalFiltered++;
                }
            }
        }

        vcb.genotypes(gc);
        writer.add(vcb.make());
    }

    @Override
    public Object onTraversalSuccess() {
        logger.warn("Total genotypes filtered: " + totalFiltered);
        return super.onTraversalSuccess();
    }

    /**
     * Closes out the new variants file.
     */
    @Override
    public void closeTool() {
        if (writer != null) {
            writer.close();
        }
    }
}
