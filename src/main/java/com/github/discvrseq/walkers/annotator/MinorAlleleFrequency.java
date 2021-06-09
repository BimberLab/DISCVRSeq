
package com.github.discvrseq.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Calculates MAF based on AF field
 *
 */
public class MinorAlleleFrequency implements InfoFieldAnnotation, StandardAnnotation {

    public static final String MAF_KEY = "MAF";

    @SuppressWarnings("unchecked")
    private List<Double> convertAF(Object obj) {
        return new ArrayList<>((List)obj);
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if ( ! vc.hasGenotypes() )
            return null;

        Map<String, Object> chrCounts = VariantContextUtils.calculateChromosomeCounts(vc, new HashMap<>(), true, Collections.emptySet());
        if (!chrCounts.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY)) {
            return null;
        }

        //this will be the frequency of each ALT allele
        List<Double> afVals;
        if (chrCounts.get(VCFConstants.ALLELE_FREQUENCY_KEY) instanceof List) {
            afVals = convertAF(chrCounts.get(VCFConstants.ALLELE_FREQUENCY_KEY));
        }
        else if (chrCounts.get(VCFConstants.ALLELE_FREQUENCY_KEY) instanceof Double) {
            afVals = new ArrayList<>();
            afVals.add((Double)chrCounts.get(VCFConstants.ALLELE_FREQUENCY_KEY));
        }
        else {
            return null;
        }

        double refAF = 1.0;
        for (Double d : afVals) {
            refAF = refAF - d;
        }
        afVals.add(refAF);

        Map<String, Object> attributeMap = new HashMap<>();
        if (afVals.size() == 1) {
            attributeMap.put(MAF_KEY, 0.0);
        }
        else {
            Collections.sort(afVals);
            attributeMap.put(MAF_KEY, afVals.get(afVals.size() - 2));
        }

        return attributeMap;
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Collections.singletonList(MAF_KEY); }

    public List<VCFCompoundHeaderLine> getDescriptions() { return Collections.singletonList(
            new VCFInfoHeaderLine(MAF_KEY, 1, VCFHeaderLineType.Float, "The minor allele frequency (frequency of second most common allele), derived from the AF field."));
    }
}
