package com.github.discvrseq.walkers.variantqc;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Molten;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;

import static com.github.discvrseq.walkers.variantqc.TableReportDescriptor.InfoFieldTableReportDescriptor.getEvalModuleSimpleName;

public class VariantEvalChild extends VariantEval {
    private final VariantQC variantQC;
    private final List<String> infoFieldEvaluators = new ArrayList<>();

    public VariantEvalChild(VariantQC variantQC, VariantQC.VariantEvalWrapper wrapper, FeatureInput<VariantContext> evals, List<String> infoFieldEvaluators){
        this.variantQC = variantQC;
        this.infoFieldEvaluators.addAll(infoFieldEvaluators);

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

    //@Override
//    protected EvaluationContext getEvaluationContext(final Set<Class<? extends VariantEvaluator>> evaluationObjects) {
//        return new ExtendedEvaluationContext(this, evaluationObjects, infoFieldEvaluators);
//    }

//    public static class ExtendedEvaluationContext extends EvaluationContext {
//        private List<String> infoFields = new ArrayList<>();
//
//        public ExtendedEvaluationContext(VariantEval walker, Set<Class<? extends VariantEvaluator>> evaluationClasses, List<String> infoFields) {
//            super(walker, evaluationClasses);
//
//            this.infoFields.addAll(infoFields);
//        }
//
//        //@Override
//        protected List<VariantEvaluator> initializeEvaluationInstances(boolean doInitialize) {
//            //TODO
//            //List<VariantEvaluator> ret = super.initializeEvaluationInstances(doInitialize);
//            List<VariantEvaluator> ret = new ArrayList<>();
//            for (String field : infoFields) {
//                ret.add(new InfoFieldEvaluator(field));
//            }
//
//            return ret;
//        }
//    }

    public static class InfoFieldEvaluator extends VariantEvaluator {
        private final String infoFieldName;

        @Molten
        Map<String, Long> counts = new HashMap<>();

        protected InfoFieldEvaluator(String infoFieldName) {
            //TODO
            //super(getEvalModuleSimpleName(infoFieldName));
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
            }
        }

        @Override
        public int getComparisonOrder() {
            return 1;
        }
    }
}
