package com.github.discvrseq.walkers.variantqc;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public class VariantEvalChild extends VariantEval {
    private final VariantQC variantQC;
    private final List<String> infoFields = new ArrayList<>();

    public VariantEvalChild(VariantQC variantQC, VariantQC.VariantEvalWrapper wrapper, FeatureInput<VariantContext> evals, List<String> infoFields){
        this.variantQC = variantQC;
        this.infoFields.addAll(infoFields);

        try
        {
            this.evals = Collections.singletonList(evals);
            this.outFile = wrapper.getOutFile();

            this.MODULES_TO_USE = new ArrayList<>(wrapper.evaluationModules);
            this.NO_STANDARD_MODULES = true;
            this.STRATIFICATIONS_TO_USE = wrapper.stratifications;
            this.NO_STANDARD_STRATIFICATIONS = true;

            //TODO: set reference??

            this.onStartup();
        }
        catch (Exception e)
        {
            throw new GATKException(e.getMessage(), e);
        }
    }

    @Override
    public List<SimpleInterval> getTraversalIntervals() {
        return variantQC.getTraversalIntervals();
    }

    @Override
    protected EvaluationContext createEvaluationContext(final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        return new ExtendedEvaluationContext(this, evaluationObjects, infoFields);
    }

    public static class ExtendedEvaluationContext extends EvaluationContext {
        private List<String> infoFields;

        public ExtendedEvaluationContext(VariantEval walker, Set<Class<? extends VariantEvaluator>> evaluationClasses, List<String> infoFields) {
            super(walker, evaluationClasses);

            this.infoFields = new ArrayList<>();
            this.infoFields.addAll(infoFields);
            for (String field : infoFields) {
                getEvaluationInstances().add(new InfoFieldEvaluator(field));
            }
        }
    }
}
