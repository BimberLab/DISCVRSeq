package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;

public class RefAlleleFrequencyArgumentCollection {
    public FeatureManager featureManager = null;

    public static interface UsesRefAlleleFrequencyArgumentCollection {
        public void setArgumentCollection(RefAlleleFrequencyArgumentCollection args);
    }

    @Argument(doc="Allele frequency source VCF", fullName = "af-source-vcf", shortName = "asv", optional = true)
    public FeatureInput<VariantContext> referenceVcf = null;

    @Argument(doc="Target INFO field key", fullName = "target-info-field-key", optional = true)
    public String targetInfoFieldKey = null;

    @Argument(doc="Source INFO field key", fullName = "source-info-field-key", optional = true)
    public String sourceInfoFieldKey = "AF";

    public void validateArguments() {
        if (referenceVcf == null) {
            throw new UserException.BadInput("Missing reference-genotypes-vcf argument");
        }

        if (targetInfoFieldKey == null) {
            throw new UserException.BadInput("Missing target-info-field-key argument");
        }
    }
}
