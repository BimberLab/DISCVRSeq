package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.PedigreeAnnotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Created by bimber on 4/20/2017.
 *
 */
public class GenotypeConcordance extends PedigreeAnnotation implements GenotypeAnnotation, GenotypeConcordanceArgumentCollection.UsesGenotypeConcordanceArgumentCollection {
    public static final String KEY = "GTD";
    public static final String D_KEY = "REF_GT";

    public GenotypeConcordanceArgumentCollection args = null;

    public GenotypeConcordance()
    {
        super((Set<String>) null);
    }

    public GenotypeConcordance(final GATKPath pedigreeFile)
    {
        super(pedigreeFile);
    }

    @Override
    public void setArgumentCollection(GenotypeConcordanceArgumentCollection args) {
        this.args = args;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(KEY, D_KEY);
    }

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (args == null) {
            throw new IllegalArgumentException("GenotypeConcordanceArgumentCollection was not set!");
        }

        if (args.referenceVcf == null) {
            throw new IllegalArgumentException("Must provide a VCF with reference genotypes!");
        }

        if (g.isFiltered() || g.isNoCall()) {
            return;
        }

        List<VariantContext> list = args.featureManager.getFeatures(args.referenceVcf, new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd()));
        if (list == null || list.isEmpty()){
            return;
        }

        for (VariantContext c : list){
            Genotype refGenotype = c.getGenotype(g.getSampleName());
            if (refGenotype != null && !refGenotype.isFiltered() && !refGenotype.isNoCall()) {
                if (!refGenotype.sameGenotype(g)) {
                    gb.attribute(KEY, "1");
                    gb.attribute(D_KEY, refGenotype.getGenotypeString());
                }
                else {
                    gb.attribute(KEY, "0");
                }

            }
        }
    }

    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFFormatHeaderLine(KEY, 1, VCFHeaderLineType.Integer, "Flags genotypes (as 1) discordant with those from the same sample/position in the provided VCF file.  Concordant genotypes are flagged a 0.  Genotypes not called in either VCF are ignored."),
                new VCFFormatHeaderLine(D_KEY, 1, VCFHeaderLineType.String, "When comparing genotypes against an alternate VCF, this will store the genotype of this sample in that alternate VCF, if discordant.  Genotypes not called in either VCF are ignored.")
        );
    }
}
