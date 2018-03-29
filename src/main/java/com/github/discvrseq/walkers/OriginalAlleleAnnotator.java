package com.github.discvrseq.walkers;

import com.github.discvrseq.annotations.OriginalAlleles;
import com.github.discvrseq.tools.DiscvrSeqDevProgramGroup;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Map;

@CommandLineProgramProperties(
        summary = "This is a fairly specialized tool that will add the current st of alleles to the INFO field as an annoation.  This was designed to preserve this information prior to liftover, during which the alleles can potentially change due to reverse complementation",
        oneLineSummary = "Annotates the current alleles",
        programGroup = DiscvrSeqDevProgramGroup.class
)
public class OriginalAlleleAnnotator extends VariantWalker {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    private VariantContextWriter writer;
    private final OriginalAlleles originalAlleles = new OriginalAlleles();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        Utils.nonNull(outFile);

        writer = createVCFWriter(new File(outFile));

        VCFHeader header = new VCFHeader(getHeaderForVariants());
        for (VCFInfoHeaderLine line : new OriginalAlleles().getDescriptions()){
            header.addMetaDataLine(line);
        }

        writer.writeHeader(header);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        Map<String, Object> annot = originalAlleles.annotate(referenceContext, variant, null);
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        if (annot != null){
            vcb.attributes(annot);
        }

        writer.add(vcb.make());
    }

    @Override
    public void closeTool(){
        writer.close();
    }
}
