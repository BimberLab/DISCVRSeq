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
import java.util.*;

public class ExtendedVariantEvalEngine extends VariantEvalEngine {

    private final List<String> infoFields = new ArrayList<>();
    private final int maxContigs;
    private final List<String> contigsToRetain;

    public ExtendedVariantEvalEngine(VariantEvalArgumentCollection variantEvalArgs, FeatureManager features, List<SimpleInterval> traversalIntervals, SAMSequenceDictionary samSequenceDictionaryForDrivingVariants, @Nullable Collection<String> samples, List<String> infoFields, int maxContigs, List<String> contigsToRetain) {
        super(variantEvalArgs, features, traversalIntervals, samSequenceDictionaryForDrivingVariants, samples);

        this.infoFields.addAll(infoFields);
        this.maxContigs = maxContigs;
        this.contigsToRetain = contigsToRetain;

        this.doValidateAndInitialize(samples);
    }

    @Override
    protected void validateAndInitialize(@Nullable Collection<String> samples) {
        // Deliberate no-op to defer initialization of EvaluationContext
    }

    private void doValidateAndInitialize(@Nullable Collection<String> samples) {
        super.validateAndInitialize(samples);
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

    @Override
    public VariantStratifier createVariantStratifier(Class<? extends VariantStratifier> clazz) {
        if (org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.Contig.class == clazz) {
            return new Contig(this, maxContigs, contigsToRetain);
        }

        return super.createVariantStratifier(clazz);
    }
}
