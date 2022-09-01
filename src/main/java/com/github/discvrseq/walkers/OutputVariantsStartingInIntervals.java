package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * As of writing, GATK SelectVariants does not support --only-output-calls-starting-in-intervals. This is a single-purpose tool that takes and input VCF, intervals, and outputs
 * only those variants that start within the supplied intervals
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar OutputVariantsStartingInIntervals \
 *     -V myFile.vcf.gz \
 *     -l 1:100-2000 \
 *     -O output.vcf.gz
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "As of writing, GATK SelectVariants does not support --only-output-calls-starting-in-intervals. This is a single-purpose tool that takes and input VCF, intervals, and outputs only those variants that start within the supplied intervals",
        oneLineSummary = "Output only variants starting within the provided intervals",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class OutputVariantsStartingInIntervals extends VariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private GATKPath outputFile;

    private VariantContextWriter vcfWriter;

    private List<SimpleInterval> intervals;

    private int totalDropped = 0;

    @Override
    public void onTraversalStart() {

        if (!hasUserSuppliedIntervals()) {
            throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required for this tool");
        }

        intervals = intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary());
        vcfWriter = createVCFWriter(outputFile);

        final VCFHeader inputVCFHeader = getHeaderForVariants();
        vcfWriter.writeHeader(inputVCFHeader);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final SimpleInterval variantStart = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getStart());
        if (intervals.stream().anyMatch(interval -> interval.contains(variantStart))) {
            vcfWriter.add(variant);
        }
        else {
            totalDropped++;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Total sites dropped: " + totalDropped);

        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
