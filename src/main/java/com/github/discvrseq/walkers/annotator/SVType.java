package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SVType implements InfoFieldAnnotation {
    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        GATKSVVCFConstants.StructuralVariantAnnotationType type = SVCallRecordUtils.inferStructuralVariantType(vc);
        if (type == null) {
            return Collections.emptyMap();
        }

        Map<String, Object> attributeMap = new HashMap<>();
        attributeMap.put(GATKSVVCFConstants.SVTYPE, type.name());

        return attributeMap;
    }

    public List<String> getKeyNames() { return Collections.singletonList(GATKSVVCFConstants.SVTYPE); }

    public List<VCFCompoundHeaderLine> getDescriptions() { return Collections.singletonList(
            GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE)
        );
    }
}
