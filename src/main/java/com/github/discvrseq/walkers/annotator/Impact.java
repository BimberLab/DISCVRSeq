
package com.github.discvrseq.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
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
public class Impact implements InfoFieldAnnotation, StandardAnnotation {

    public static final String IMPACT_KEY = "IMPACT";

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (!vc.hasAttribute("ANN")) {
            return null;
        }

        List<String> anns = vc.getAttributeAsStringList("ANN", null);
        if (anns == null) {
            return null;
        }

        Set<String> impacts = new HashSet<>();
        for (String ann : anns) {
            String[] split = ann.split("\\|");
            if (split.length > 7 && "protein_coding".equals(split[7])) {
                impacts.add(split[2]);
            }
        }

        Map<String, Object> attributeMap = new HashMap<>();
        for (String val : Arrays.asList("HIGH", "LOW", "MODERATE")) {
            if (impacts.contains(val)) {
                attributeMap.put(IMPACT_KEY, val);
                break;
            }
        }

        if (attributeMap.isEmpty()) {
            return null;
        }

        return attributeMap;
    }

    public List<String> getKeyNames() { return Collections.singletonList(IMPACT_KEY); }

    public List<VCFCompoundHeaderLine> getDescriptions() { return Collections.singletonList(
            new VCFInfoHeaderLine(IMPACT_KEY, 1, VCFHeaderLineType.String, "The highest impact annotation provided by SnpEff, limited to protein_coding features"));
    }
}
