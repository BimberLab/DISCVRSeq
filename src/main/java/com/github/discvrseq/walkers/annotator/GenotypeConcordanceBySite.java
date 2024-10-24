package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.PedigreeAnnotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by bimber on 4/20/2017.
 *
 */
public class GenotypeConcordanceBySite extends PedigreeAnnotation implements InfoFieldAnnotation, GenotypeConcordanceArgumentCollection.UsesGenotypeConcordanceArgumentCollection {
    public static final String DISCORD_KEY = "GTD";
    public static final String CONCORD_KEY = "GTC";

    public GenotypeConcordanceArgumentCollection args = null;

    public GenotypeConcordanceBySite()
    {
        super((Set<String>) null);
    }

    public GenotypeConcordanceBySite(final GATKPath pedigreeFile)
    {
        super(pedigreeFile);
    }

    @Override
    public void setArgumentCollection(GenotypeConcordanceArgumentCollection args) {
        this.args = args;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(DISCORD_KEY, CONCORD_KEY);
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (args == null) {
            throw new IllegalArgumentException("GenotypeConcordanceArgumentCollection was not set!");
        }

        if (args.referenceVcf == null) {
            throw new IllegalArgumentException("Must provide a VCF with reference genotypes!");
        }

        List<VariantContext> referenceVCs = args.featureManager.getFeatures(args.referenceVcf, new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd())).stream().filter(refVC -> refVC.getStart() == vc.getStart()).toList();
        if (referenceVCs.isEmpty()){
            return null;
        }

        if (referenceVCs.size() > 1) {
            throw new IllegalArgumentException("More than one reference found for site: " + vc.getContig() + "-" + vc.getStart());
        }

        VariantContext refVC  = referenceVCs.get(0);

        AtomicInteger discord = new AtomicInteger(0);
        AtomicInteger concord = new AtomicInteger(0);
        Iterator<Genotype> it = vc.getGenotypes().iterator();
        it.forEachRemaining(g -> {
            if (!g.isFiltered() && !g.isNoCall()) {
                Genotype refGenotype = refVC.getGenotype(g.getSampleName());
                if (refGenotype != null && !refGenotype.isFiltered() && !refGenotype.isNoCall()) {
                    if (!refGenotype.sameGenotype(g)) {
                        discord.getAndIncrement();
                    }
                    else {
                        concord.getAndIncrement();
                    }
                }
            }
        });

        Map<String,Object> attributeMap = new HashMap<>(2);
        attributeMap.put(DISCORD_KEY, discord.get());
        attributeMap.put(CONCORD_KEY, concord.get());

        return attributeMap;
    }

    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFInfoHeaderLine(DISCORD_KEY, 1, VCFHeaderLineType.Integer, "The total number of genotypes discordant with those from the same sample/position in the provided VCF file.  Genotypes not called in either VCF are ignored."),
                new VCFInfoHeaderLine(CONCORD_KEY, 1, VCFHeaderLineType.Integer, "The total number of genotypes concordant with those from the same sample/position in the provided VCF file.  Genotypes not called in either VCF are ignored.")
        );
    }
}
