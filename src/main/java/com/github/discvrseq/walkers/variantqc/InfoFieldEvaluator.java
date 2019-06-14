package com.github.discvrseq.walkers.variantqc;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Molten;

import java.util.HashMap;
import java.util.Map;

import static com.github.discvrseq.walkers.variantqc.TableReportDescriptor.InfoFieldTableReportDescriptor.getEvalModuleSimpleName;

@Analysis(
        name = "Info Field Summary",
        description = "Summary of the specified INFO field",
        molten = true
)
public class InfoFieldEvaluator extends VariantEvaluator {
    private final String infoFieldName;

    @Molten (
        variableName = "Value",
        valueName = "Count"
    )
    public Map<String, Long> counts = new HashMap<>();
    private long total = 0;

    protected InfoFieldEvaluator(String infoFieldName) {
        super(getEvalModuleSimpleName(infoFieldName));
        this.infoFieldName = infoFieldName;
    }

    @Override
    public void update1(VariantContext eval, ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext) {
        if (eval != null && eval.hasAttribute(infoFieldName)) {
            Object val = eval.getAttribute(infoFieldName);
            if (val != null) {
                String stringVal = val.toString();

                Long count = counts.getOrDefault(stringVal, 0L);
                count += 1;

                counts.put(stringVal, count);
            }

            total += 1;
        }
    }

    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void finalizeEvaluation() {
        if (counts.isEmpty()) {
            counts.put("Empty/Blank", total);
        }
    }
}
