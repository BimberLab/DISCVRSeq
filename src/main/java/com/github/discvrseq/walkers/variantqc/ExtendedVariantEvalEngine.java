package com.github.discvrseq.walkers.variantqc;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

public class ExtendedVariantEvalEngine extends VariantEvalEngine {

    private final List<String> infoFields = new ArrayList<>();
    private final int maxContigs;
    private final List<String> contigsToRetain;

    private ExtendedVariantEvalEngine(VariantEvalArgumentCollection variantEvalArgs, FeatureManager features, List<SimpleInterval> traversalIntervals, SAMSequenceDictionary samSequenceDictionaryForDrivingVariants, @Nullable Collection<String> samples, List<String> infoFields, int maxContigs, List<String> contigsToRetain) {
        super(variantEvalArgs, features, traversalIntervals, samSequenceDictionaryForDrivingVariants, samples);

        this.infoFields.addAll(infoFields);
        this.maxContigs = maxContigs;
        this.contigsToRetain = contigsToRetain;
    }

    public static ExtendedVariantEvalEngine create(VariantEvalArgumentCollection variantEvalArgs, FeatureManager features, List<SimpleInterval> traversalIntervals, SAMSequenceDictionary samSequenceDictionaryForDrivingVariants, @Nullable Collection<String> samples, List<String> infoFields, int maxContigs, List<String> contigsToRetain)
    {
        ExtendedVariantEvalEngine ee = new ExtendedVariantEvalEngine(variantEvalArgs, features, traversalIntervals, samSequenceDictionaryForDrivingVariants, samples, infoFields, maxContigs, contigsToRetain);
        ee.doValidateAndInitializeAfterConstruction(samples);

        return ee;
    }

    @Override
    protected void validateAndInitialize(@Nullable Collection<String> samples) {
        // Deliberate no-op to defer initialization of EvaluationContext
    }

    private void doValidateAndInitializeAfterConstruction(@Nullable Collection<String> samples) {
        super.validateAndInitialize(samples);
    }

    @Override
    protected EvaluationContext createEvaluationContext(final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
        return ExtendedEvaluationContext.create(this, evaluationObjects, infoFields);
    }

    public static class ExtendedEvaluationContext extends EvaluationContext {
        private ExtendedEvaluationContext(VariantEvalEngine engine, Set<Class<? extends VariantEvaluator>> evaluationClasses, List<String> infoFields) {
            super(engine, evaluationClasses);
        }

        public static ExtendedEvaluationContext create(VariantEvalEngine engine, Set<Class<? extends VariantEvaluator>> evaluationClasses, List<String> infoFields)
        {
            ExtendedEvaluationContext eec = new ExtendedEvaluationContext(engine, evaluationClasses, infoFields);
            eec.initializeEvaluationInstances(engine, infoFields);

            return eec;
        }

        // Separate to avoid this-escape warning
        protected void initializeEvaluationInstances(VariantEvalEngine engine, List<String> infoFields)
        {
            for (String field : infoFields) {
                getEvaluationInstances().add(new InfoFieldEvaluator(engine, field));
            }
        }
    }

    @Override
    public VariantStratifier createVariantStratifier(Class<? extends VariantStratifier> clazz) {
        if (org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.Contig.class == clazz) {
            return new Contig(this, maxContigs, contigsToRetain);
        }

        return super.createVariantStratifier(clazz);
    }
}
