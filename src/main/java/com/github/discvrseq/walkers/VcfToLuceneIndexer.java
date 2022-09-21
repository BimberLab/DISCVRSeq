package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.github.discvrseq.walkers.annotator.DiscvrVariantAnnotator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;

import java.util.*;

/**
 * This tool will accept a VCF, and iterate the variants and write the results to a lucene search index.
 * It accepts a list of INFO annotations that will be indexed per site.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar VcfToLuceneIndexer
 *     -V myVcf.gz \
 *     -O /directory/for/output/
 *     -IF AF
 *     -AN SampleList
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool will accept a VCF, and iterate the variants and write the results to a lucene search index",
        oneLineSummary = "Index a VCF using Apache Lucene",
        programGroup = DiscvrSeqProgramGroup.class
)
public class VcfToLuceneIndexer extends VariantWalker {
    @Argument(doc="Folder where the index will be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath outDir = null;

    @Argument(doc="Info fields to index", fullName = "info-field", shortName = "IF", optional = true)
    List<String> infoFieldsToIndex;

    @Argument(doc="Annotations to apply", fullName = "annotation-name", shortName = "AN", optional = true)
    List<String> annotationClassNames;

    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        // This will eventually be used if we need to implement custom InfoFieldAnnotation classes
        return Collections.singletonList(new DiscvrVariantAnnotator.DiscvrAnnotationPluginDescriptor());
    }

    @Override
    public boolean useVariantAnnotations() {
        return true;
    }

    private final List<InfoFieldAnnotation> annotationsToUse = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        //IOUtil.assertDirectoryIsWritable(outDir.toPath());

        // This is where you'd sanity check the set of annotations is valid, etc.
        // Note: we can implement custom InfoFieldAnnotation classes for operations like indexing all

        // Also initialize whatever kind of lucene writer is needed.

        // NOTE: as far as listing variable samples, see the GATK SampleList annotation
        DiscvrVariantAnnotator.DiscvrAnnotationPluginDescriptor plugin = getCommandLineParser().getPluginDescriptor(DiscvrVariantAnnotator.DiscvrAnnotationPluginDescriptor.class);
        for (String className : annotationClassNames) {
            Class<?> clazz = plugin.getClassForPluginHelp(className);

            try {
                Annotation a = new DiscvrVariantAnnotator.DiscvrAnnotationPluginDescriptor().createInstanceForPlugin(clazz);
                if (!(a instanceof InfoFieldAnnotation)) {
                    throw new GATKException("All custom annotations must be InfoFieldAnnotation: " + className);
                }
                annotationsToUse.add((InfoFieldAnnotation) a);
            }
            catch (InstantiationException | IllegalAccessException e) {
                throw new GATKException("Unable to create annotation: " + className);
            }
        }
    }

    private long sites = 0;

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        //TODO: per variant, build the document and write to the index

        Map<String, Object> toIndex = new HashMap<>();

        // Standard fields also included:
        toIndex.put("contig", variant.getContig());
        toIndex.put("start", variant.getStart());
        toIndex.put("end", variant.getEnd());
        // TODO: maybe others?


        for (String infoField : infoFieldsToIndex) {
            if (variant.hasAttribute(infoField) && variant.getAttribute(infoField) != null) {
                toIndex.put(infoField, variant.getAttribute(infoField));
            }
        }

        // TODO: maybe warn/throw if there are conflicting keys?
        for (InfoFieldAnnotation a : annotationsToUse) {
            toIndex.putAll(a.annotate(referenceContext, variant, null));
        }

        sites++;

        //TODO: actually index that
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Indexing complete, total sites indexed: " + sites);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        //TODO: if you need to close the Writer or anything, do that here
    }
}
