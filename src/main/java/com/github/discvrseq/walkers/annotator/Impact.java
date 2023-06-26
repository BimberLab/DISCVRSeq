
package com.github.discvrseq.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
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
    public static final String IMPACT_GENES_KEY = "HIG";

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (!vc.hasAttribute("ANN")) {
            return null;
        }

        List<String> anns = vc.getAttributeAsStringList("ANN", null);
        if (anns == null) {
            return null;
        }

        Map<String, Object> attributeMap = new HashMap<>();
        Map<String, List<String[]>> annByAllele = new HashMap<>();
        for (String ann : anns) {
            String[] split = ann.split("\\|");
            if (!annByAllele.containsKey(split[0])) {
                annByAllele.put(split[0], new ArrayList<>());
            }

            annByAllele.get(split[0]).add(split);
        }

        Map<Allele, String> impactMap = new HashMap<>();
        Map<Allele, String> impactGeneMap = new HashMap<>();

        for (Allele a : vc.getAlternateAlleles()) {
            if (!annByAllele.containsKey(a.getBaseString())) {
                continue;
            }

            Set<String> impacts = new HashSet<>();
            Set<String> hig = new HashSet<>();
            for (String[] split : annByAllele.get(a.getBaseString())) {
                if (split.length > 7 && "protein_coding".equals(split[7])) {
                    impacts.add(split[2]);

                    if ("HIGH".equals(split[2])) {
                        hig.add(split[3]);
                    }
                }
            }

            for (String val : Arrays.asList("HIGH", "LOW", "MODERATE")) {
                if (impacts.contains(val)) {
                    impactMap.put(a, val);
                    break;
                }
            }

            if (!hig.isEmpty()) {
                impactGeneMap.put(a, StringUtils.join(hig, "|"));
            }
        }

        if (!impactMap.isEmpty()) {
            attributeMap.put(IMPACT_KEY, vc.getAlternateAlleles().stream().map(a -> impactMap.getOrDefault(a, "")).toList());
        }

        if (!impactMap.isEmpty()) {
            attributeMap.put(IMPACT_GENES_KEY, vc.getAlternateAlleles().stream().map(a -> impactGeneMap.getOrDefault(a, "")).toList());
        }

        if (attributeMap.isEmpty()) {
            return null;
        }

        return attributeMap;
    }

    public List<String> getKeyNames() { return Arrays.asList(IMPACT_KEY, IMPACT_GENES_KEY); }

    public List<VCFCompoundHeaderLine> getDescriptions() { return Arrays.asList(
            new VCFInfoHeaderLine(IMPACT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "The highest impact annotation provided by SnpEff, limited to protein_coding features"),
            new VCFInfoHeaderLine(IMPACT_GENES_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "A comma-separated list of any overlapping genes with high-impact effects on protein coding, identified by SnpEff")
        );
    }
}
