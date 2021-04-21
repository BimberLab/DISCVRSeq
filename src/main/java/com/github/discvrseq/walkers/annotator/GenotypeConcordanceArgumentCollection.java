package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;

public class GenotypeConcordanceArgumentCollection {
    @Argument(doc="Reference genotypes VCF", fullName = "reference-genotypes-vcf", shortName = "rg", optional = true)
    public FeatureInput<VariantContext> referenceVcf = null;

    public static interface UsesGenotypeConcordanceArgumentCollection {
        public void setArgumentCollection(GenotypeConcordanceArgumentCollection args);
    }

    public void validateArguments() {
        if (referenceVcf == null) {
            throw new UserException.BadInput("Missing reference-genotypes-vcf argument");
        }
    }

    public FeatureManager featureManager = null;
}
