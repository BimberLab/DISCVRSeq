package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Created by bimber on 4/20/2017.
 *
 */
public class RefAlleleFrequency implements InfoFieldAnnotation, RefAlleleFrequencyArgumentCollection.UsesRefAlleleFrequencyArgumentCollection {

    public RefAlleleFrequencyArgumentCollection args = null;

    public RefAlleleFrequency() {

    }

    @Override
    public void setArgumentCollection(RefAlleleFrequencyArgumentCollection args) {
        this.args = args;
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(args.targetInfoFieldKey);
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (args == null) {
            throw new IllegalArgumentException("RefAlleleFrequencyArgumentCollection was not set!");
        }

        List<VariantContext> list = args.featureManager.getFeatures(args.referenceVcf, new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd()));
        if (list == null || list.isEmpty()){
            return null;
        }

        // Require the ref allele to be identical:
        list = list.stream().filter(x -> vc.getReference().equals(x.getReference())).filter(x -> x.hasAttribute(args.sourceInfoFieldKey) && x.getAttribute(args.sourceInfoFieldKey) != null).toList();
        if (list.isEmpty()){
            return null;
        }
        else if (list.size() > 1) {
            // TODO: consider logging this?
            return null;
        }

        VariantContext vc2 = list.get(0);
        List<Double> afs = vc2.getAttributeAsDoubleList(args.sourceInfoFieldKey, 0.0);
        if (afs.size() != vc2.getAlternateAlleles().size()) {
            throw new UserException.BadInput("Source field and alt alleles are not the same size: " + vc2.toStringWithoutGenotypes());
        }

        List<Double> toAdd = vc.getAlternateAlleles().stream().map(a -> {
            if (vc2.hasAllele(a)) {
                return afs.get(vc2.getAlleleIndex(a) - 1);
            }
            else {
                return 0.0;
            }
        }).toList();

        Map<String, Object> ret = new HashMap<>();
        ret.put(args.targetInfoFieldKey, toAdd);

        return ret;
    }

    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return Collections.singletonList(
                new VCFInfoHeaderLine(args.targetInfoFieldKey, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "")
        );
    }
}
