package com.github.discvrseq.walkers.variantqc;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import javax.annotation.Nullable;
import java.util.*;

public class ExtendedVariantEvalEngine extends VariantEvalEngine {

    private final List<String> infoFields = new ArrayList<>();

    public ExtendedVariantEvalEngine(VariantEvalArgumentCollection variantEvalArgs, FeatureManager features, List<SimpleInterval> traversalIntervals, SAMSequenceDictionary samSequenceDictionaryForDrivingVariants, @Nullable Collection<String> samples, List<String> infoFields) {
        super(variantEvalArgs, features, traversalIntervals, samSequenceDictionaryForDrivingVariants, samples, true);

        this.infoFields.addAll(infoFields);
        validateAndInitialize(samples);
    }

    @Override
    protected EvaluationContext createEvaluationContext(final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        return new ExtendedEvaluationContext(this, evaluationObjects, infoFields);
    }

    public static class ExtendedEvaluationContext extends EvaluationContext {
        public ExtendedEvaluationContext(VariantEvalEngine engine, Set<Class<? extends VariantEvaluator>> evaluationClasses, List<String> infoFields) {
            super(engine, evaluationClasses);

            for (String field : infoFields) {
                getEvaluationInstances().add(new InfoFieldEvaluator(engine, field));
            }
        }
    }
}
