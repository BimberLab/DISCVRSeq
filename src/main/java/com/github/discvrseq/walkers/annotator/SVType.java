package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

public class SVType implements InfoFieldAnnotation {
    private static final Set<String> VALID_TYPES = new HashSet<>(Arrays.asList(GATKSVVCFConstants.StructuralVariantAnnotationType.values()).stream()
            .map(GATKSVVCFConstants.StructuralVariantAnnotationType::name).collect(Collectors.toList()));

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (!vc.isVariant()) {
            return Collections.emptyMap();
        }

        if (vc.getType() == VariantContext.Type.SNP || vc.getType() == VariantContext.Type.MNP) {
            return Collections.emptyMap();
        }

        Set<GATKSVVCFConstants.StructuralVariantAnnotationType> types = new HashSet<>();
        for (Allele a : vc.getAlternateAlleles()) {
            if (vc.isSimpleInsertion()) {
                types.add(GATKSVVCFConstants.StructuralVariantAnnotationType.INS);
            }
            else if (vc.isSimpleDeletion()) {
                types.add(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
            }
            if (a.isSymbolic()) {
                final String alleleType = a.getDisplayString().replace("<", "").replace(">", "");
                if (VALID_TYPES.contains(alleleType)) {
                    types.add(GATKSVVCFConstants.StructuralVariantAnnotationType.valueOf(alleleType));
                } else {
                    throw new IllegalArgumentException("Could not find a valid SV type for variant " + vc.getID());
                }
            }
            else if (a.isBreakpoint()){
                types.add(GATKSVVCFConstants.StructuralVariantAnnotationType.BND);
            }
            else if (a.length() > vc.getReference().length()) {
                types.add(GATKSVVCFConstants.StructuralVariantAnnotationType.INS);
            }
            else if (a.length() < vc.getReference().length()) {
                types.add(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
            }
        }

        Map<String, Object> attributeMap = new HashMap<>();
        if (types.size() == 1) {
            String type = types.stream().map(GATKSVVCFConstants.StructuralVariantAnnotationType::name).findFirst().get();
            attributeMap.put(GATKSVVCFConstants.SVTYPE, type);
        }
        else if (types.size() > 1) {
            attributeMap.put(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CPX);
        }

        return attributeMap;
    }

    public List<String> getKeyNames() { return Collections.singletonList(GATKSVVCFConstants.SVTYPE); }

    public List<VCFCompoundHeaderLine> getDescriptions() { return Collections.singletonList(
            GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE)
        );
    }
}
