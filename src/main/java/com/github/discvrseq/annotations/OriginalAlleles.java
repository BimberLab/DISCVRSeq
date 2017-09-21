package com.github.discvrseq.annotations;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * This annotation is primarily intended to be used for situations where a variant
 * is lifted to another genome, which could result in reverse complementation or shifting of alleles
 */
public class OriginalAlleles extends InfoFieldAnnotation {
    public static final String KEY = "OA";

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, ReadLikelihoods<Allele> likelihoods) {
        List<String> alleles = vc.getAlleles().stream().map(a -> a.getBaseString()).collect(Collectors.toList());
        return Collections.singletonMap(KEY, StringUtils.join(alleles, ","));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(new VCFInfoHeaderLine(KEY, VCFHeaderLineCount.R, VCFHeaderLineType.String, "This contains a list of the original alleles at this site.  It is primarily intended to store information prior to liftover."));
    }
}
