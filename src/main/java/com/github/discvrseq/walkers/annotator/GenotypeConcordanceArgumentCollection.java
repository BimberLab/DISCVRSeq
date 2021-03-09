package com.github.discvrseq.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;

public class GenotypeConcordanceArgumentCollection {
    @Argument(doc="Reference genotypes VCF", fullName = "reference-genotypes-vcf", shortName = "rg", optional = true)
    public FeatureInput<VariantContext> referenceVcf = null;
}
