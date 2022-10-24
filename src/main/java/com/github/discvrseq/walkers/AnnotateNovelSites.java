package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * This is a fairly specialized tool, designed to take one reference VCF and use this to annotate another VCF to report which sites are not present in the original VCF
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar AnnotateNovelSites \
 *     -V input.vcf.gz \
 *     -rv referenceSites.vcf \
 *     -ns novelSitesOutput.vcf \
 *     -ms missingSitesOutput.vcf \
 *     -an mGAPV \
 *     -ad 'The first mGAP version where variants at this site appeared' \
 *     -av v2.0 \
 *     -O output.annotated.vcf.gz
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This is a fairly specialized tool, designed to take one reference VCF and use this to annotate another VCF to report which sites are not present in the original VCF",
        oneLineSummary = "Annotates a VCF based on whether sites are present in a reference VCF",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class AnnotateNovelSites extends ExtendedMultiVariantWalkerGroupedOnStart {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath outFile = null;

    @Argument(doc="Novel Sites VCF", fullName = "novel-sites", shortName = "ns", optional = true)
    public GATKPath novelSitesOutput = null;

    @Argument(doc="Missing Sites VCF", fullName = "missing-sites", shortName = "ms", optional = true)
    public GATKPath droppedSitesOutput = null;

    @Argument(doc="Novel Site Annotation Value", fullName = "novel-site-annotation", shortName = "av", optional = false)
    public String novelSiteAnnotationValue = null;

    @Argument(doc="Existing Site Annotation Default Value", fullName = "existing-site-annotation-default", shortName = "dv", optional = true)
    public String existingSiteDefaultAnnotationValue = null;

    @Argument(doc="Novel Site Annotation Name", fullName = "novel-site-annotation-name", shortName = "an", optional = false)
    public String novelSiteAnnotationName = null;

    @Argument(doc="Novel Site Annotation Description", fullName = "novel-site-annotation-description", shortName = "ad", optional = false)
    public String novelSiteAnnotationDescription = null;

    private VariantContextWriter writer;
    private VariantContextWriter novelSitesWriter = null;
    private VariantContextWriter droppedSitesWriter = null;

    private long novelSites = 0L;
    private long droppedSites = 0L;

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        if (getMultiVariantInputArgumentCollection().getDrivingVariantPaths().size() > 1) {
            throw new IllegalArgumentException("Cannot provide more than one input VCF (-V)");
        }

        writer = createVCFWriter(outFile.toPath());

        VCFHeader sourceHeader = (VCFHeader)getHeaderForFeatures(getDrivingVariantsFeatureInputs().get(0));

        VCFHeader header = new VCFHeader(sourceHeader.getMetaDataInInputOrder(), sourceHeader.getSampleNamesInOrder());

        if (!header.hasInfoLine(novelSiteAnnotationName)) {
            header.addMetaDataLine(new VCFInfoHeaderLine(novelSiteAnnotationName, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, novelSiteAnnotationDescription));
        }

        writer.writeHeader(header);

        VCFHeader headerSamples = new VCFHeader(header.getMetaDataInInputOrder());
        if (novelSitesOutput != null) {
            novelSitesWriter = createVCFWriter(novelSitesOutput.toPath());
            novelSitesWriter.writeHeader(headerSamples);
        }

        if (droppedSitesOutput != null) {
            droppedSitesWriter = createVCFWriter(droppedSitesOutput.toPath());
            droppedSitesWriter.writeHeader(headerSamples);
        }
    }

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        Map<FeatureInput<VariantContext>, List<VariantContext>> variants = groupVariantsByFeatureInput(variantContexts);
        List<VariantContext> refVariants = variants.get(getAnnotateNovelSitesArgumentCollection().referenceVcf);

        List<VariantContext> inputVariants = variants.get(getDrivingVariantsFeatureInputs().get(0));
        if (inputVariants == null || inputVariants.isEmpty()) {
            if (droppedSitesWriter != null) {
                droppedSites++;

                refVariants.forEach(vc -> {
                    VariantContextBuilder vcb = new VariantContextBuilder(vc).genotypes();
                    droppedSitesWriter.add(vcb.make());
                });
            }

            return;
        }

        for (VariantContext vc : inputVariants) {
            VariantContextBuilder vcb = new VariantContextBuilder(vc);

            // Site is novel:
            if (refVariants == null || refVariants.isEmpty()) {
                vcb.attribute(novelSiteAnnotationName, novelSiteAnnotationValue);
                writer.add(vcb.make());

                if (novelSitesWriter != null) {
                    novelSites++;
                    novelSitesWriter.add(vcb.genotypes().make());
                }
            }
            else {
                boolean foundMatchingRef = false;
                boolean isNovel = false;
                for (VariantContext refVariant : refVariants) {
                    // Prior site exists, alleles identical. No action needed:
                    if (refVariant.getReference().equals(vc.getReference())) {
                        if (refVariant.getAlternateAlleles().equals(vc.getAlternateAlleles())) {
                            List<String> val = refVariant.getAttributeAsStringList(novelSiteAnnotationName, existingSiteDefaultAnnotationValue);
                            if (val != null && !val.isEmpty()) {
                                vcb.attribute(novelSiteAnnotationName, new ArrayList<>(val));
                            }
                        }
                        else {
                            // Alleles are different, so annotate with both
                            Set<String> anns = new LinkedHashSet<>(refVariant.getAttributeAsStringList(novelSiteAnnotationName, existingSiteDefaultAnnotationValue));
                            anns.add(novelSiteAnnotationValue);
                            isNovel = true;

                            vcb.attribute(novelSiteAnnotationName, new ArrayList<>(anns));
                        }
                        foundMatchingRef = true;
                        break;
                    }
                }

                if (!foundMatchingRef) {
                    vcb.attribute(novelSiteAnnotationName, novelSiteAnnotationValue);
                    isNovel = true;
                }

                vc = vcb.make();
                writer.add(vc);

                if (isNovel && novelSitesWriter != null) {
                    novelSites++;
                    novelSitesWriter.add(vcb.genotypes().make());
                }
            }
        }
    }

    @Override
    public void closeTool(){
        writer.close();
        if (novelSitesWriter != null) {
            novelSitesWriter.close();
        }

        if (droppedSitesWriter != null) {
            droppedSitesWriter.close();
        }

        logger.info("Total novel sites: " + novelSites);
        logger.info("Total missing/dropped sites: " + droppedSites);
    }

    @Override
    protected AnnotateNovelSitesArgumentCollection getMultiVariantInputArgumentCollection() {
        return new AnnotateNovelSitesArgumentCollection();
    }

    private AnnotateNovelSitesArgumentCollection getAnnotateNovelSitesArgumentCollection() {
        return (AnnotateNovelSitesArgumentCollection)multiVariantInputArgumentCollection;
    }

    private static final class AnnotateNovelSitesArgumentCollection extends MultiVariantInputArgumentCollection.DefaultMultiVariantInputArgumentCollection {
        private static final long serialVersionUID = 1L;

        @Argument(doc="Reference VCF", fullName = "ref-vcf", shortName = "rv", optional = false)
        public FeatureInput<VariantContext> referenceVcf = null;

        @Override
        public List<GATKPath> getDrivingVariantPaths() {
            List<GATKPath> ret = new ArrayList<>(super.getDrivingVariantPaths());
            ret.add(referenceVcf);

            return ret;
        }
    }
}
