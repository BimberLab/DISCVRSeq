package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.PedigreeAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.samples.SampleDB;

import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Created by bimber on 3/13/2017.
 *
 */
public class MendelianViolationBySample extends PedigreeAnnotation implements GenotypeAnnotation, MendelianViolationArgumentCollection.UsesMendelianViolationArgumentCollection {
    private SampleDB sampleDB = null;

    public static final String MV_KEY = "MV";

    public MendelianViolationArgumentCollection args = null;

    public MendelianViolationBySample()
    {
        super((Set<String>) null);
    }

    public MendelianViolationBySample(final GATKPath pedigreeFile){
        super(pedigreeFile);
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(MV_KEY);
    }

    @Override
    public void setArgumentCollection(MendelianViolationArgumentCollection args) {
        this.args = args;
    }

    @Override
    public void annotate(ReferenceContext referenceContext, VariantContext variantContext, Genotype genotype, GenotypeBuilder gb, AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods) {
        for (String key : getKeyNames())
        {
            gb.attribute(key, null);
        }

        if (sampleDB != null) {
            int totalViolations = MendelianViolationCount.countViolations(sampleDB.getSample(genotype.getSampleName()), variantContext, args.minGenotypeQuality);
            gb.attribute(MV_KEY, totalViolations);
        }
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(new VCFFormatHeaderLine(MV_KEY, 1, VCFHeaderLineType.Integer, "Number of mendelian violations observed for this sample."));
    }
}
