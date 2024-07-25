package com.github.discvrseq.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;

import java.util.*;

abstract public class ExtendedMultiVariantWalkerGroupedOnStart extends MultiVariantWalkerGroupedOnStart {
    // maintain the mapping of source name (from VC) to FeatureInput name
    private Map<String, FeatureInput<VariantContext>> drivingVariantSourceMap = null;

    protected Map<String, FeatureInput<VariantContext>> getDrivingVariantSourceMap() {
        if (drivingVariantSourceMap == null) {
            // Cache map of source name -> FeatureInput
            drivingVariantSourceMap = new HashMap<>();
            getDrivingVariantsFeatureInputs().forEach(x -> drivingVariantSourceMap.put(x.getName(), x));
        }

        return drivingVariantSourceMap;
    }

    protected Map<FeatureInput<VariantContext>, List<VariantContext>> groupVariantsByFeatureInput(final List<VariantContext> variants) {
        final Map<FeatureInput<VariantContext>, List<VariantContext>> byFeatureInput = new HashMap<>();
        variants.forEach(vc -> byFeatureInput.compute(getDrivingVariantSourceMap().get(vc.getSource()),
                (k, v) -> {
                    final List<VariantContext> variantList = v == null ? new ArrayList<>() : v;
                    variantList.add(vc);
                    return variantList;
                }
        ));

        for (FeatureInput<VariantContext> fi : getDrivingVariantsFeatureInputs()) {
            if (!byFeatureInput.containsKey(fi)) {
                byFeatureInput.put(fi, Collections.emptyList());
            }
            else {
                byFeatureInput.get(fi).sort(new VariantContextComparator(getBestAvailableSequenceDictionary()));
            }
        }

        return byFeatureInput;
    }
}
