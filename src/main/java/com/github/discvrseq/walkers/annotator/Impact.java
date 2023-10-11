
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
    public static final String OVERLAPPING_GENES_KEY = "OG";
    public static final String EFFECT_KEY = "VE";

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
        Map<Allele, String> allGeneMap = new HashMap<>();
        Map<Allele, String> effectMap = new HashMap<>();

        for (Allele a : vc.getAlternateAlleles()) {
            // Note: use displayString to handle symbolic alleles
            if (!annByAllele.containsKey(a.getDisplayString())) {
                continue;
            }

            Set<String> impacts = new HashSet<>();
            Set<String> hig = new HashSet<>();
            Set<String> og = new HashSet<>();
            Set<String> effects = new HashSet<>();
            for (String[] split : annByAllele.get(a.getDisplayString())) {
                if (split.length > 7 && "protein_coding".equals(split[7])) {
                    impacts.add(split[2]);

                    if (!split[3].isEmpty()) {
                        og.add(split[3]);
                        if ("HIGH".equals(split[2])) {
                            hig.add(split[3]);
                        }
                    }
                }

                if (!split[1].isEmpty()) {
                    effects.add(split[1]);
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

            if (!og.isEmpty()) {
                allGeneMap.put(a, StringUtils.join(og, "|"));
            }

            if (!effects.isEmpty()) {
                effectMap.put(a, StringUtils.join(effects, "|"));
            }
        }

        if (!impactMap.isEmpty()) {
            attributeMap.put(IMPACT_KEY, vc.getAlternateAlleles().stream().map(a -> impactMap.getOrDefault(a, "")).toList());
        }

        if (!impactGeneMap.isEmpty()) {
            attributeMap.put(IMPACT_GENES_KEY, vc.getAlternateAlleles().stream().map(a -> impactGeneMap.getOrDefault(a, "")).toList());
        }

        if (!allGeneMap.isEmpty()) {
            attributeMap.put(OVERLAPPING_GENES_KEY, vc.getAlternateAlleles().stream().map(a -> allGeneMap.getOrDefault(a, "")).toList());
        }

        if (!effectMap.isEmpty()) {
            attributeMap.put(EFFECT_KEY, vc.getAlternateAlleles().stream().map(a -> effectMap.getOrDefault(a, "")).toList());
        }

        if (attributeMap.isEmpty()) {
            return null;
        }

        return attributeMap;
    }

    public List<String> getKeyNames() { return Arrays.asList(IMPACT_KEY, IMPACT_GENES_KEY, OVERLAPPING_GENES_KEY, EFFECT_KEY); }

    public List<VCFCompoundHeaderLine> getDescriptions() { return Arrays.asList(
            new VCFInfoHeaderLine(IMPACT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "The highest impact annotation provided by SnpEff, limited to protein_coding features"),
            new VCFInfoHeaderLine(IMPACT_GENES_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "A comma-separated list of any overlapping genes with high-impact effects on protein coding, identified by SnpEff"),
            new VCFInfoHeaderLine(OVERLAPPING_GENES_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "A comma-separated list of any overlapping genes, identified by SnpEff. This includes variants immediately upstream/downstream of the coding region."),
            new VCFInfoHeaderLine(EFFECT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.String, "A comma-separated list of predicted variant effects, generated by SnpEff")
        );
    }
}
