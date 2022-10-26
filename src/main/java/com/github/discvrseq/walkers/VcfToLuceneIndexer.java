package com.github.discvrseq.walkers;

import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.store.Directory;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.FloatPoint;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.TextField;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.github.discvrseq.walkers.annotator.DiscvrVariantAnnotator;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * This tool accepts a VCF, iterates the variants and writes the results to a lucene search index.
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
        summary = "This tool accepts a VCF, iterates the variants, and writes the results to a lucene search index",
        oneLineSummary = "Index a VCF using Apache Lucene",
        programGroup = DiscvrSeqProgramGroup.class
)
public class VcfToLuceneIndexer extends VariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outDir = null;

    @Argument(doc="Info fields to index", fullName = "info-field", shortName = "IF", optional = true)
    List<String> infoFieldsToIndex;

    @Argument(doc="Annotations to apply", fullName = "annotation-name", shortName = "AN", optional = true)
    List<String> annotationClassNames;

    IndexWriter writer;

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
        StandardAnalyzer analyzer = new StandardAnalyzer();

        Directory index;
        try {
            index = FSDirectory.open(outDir.toPath());
        } catch(IOException e) {
            analyzer.close();
            throw new GATKException(e.getMessage(), e);
        }

        IndexWriterConfig config = new IndexWriterConfig(analyzer);

        try {
            writer = new IndexWriter(index, config);
        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        // This is where you'd sanity check the set of annotations is valid, etc.
        // Note: we can implement custom InfoFieldAnnotation classes for operations like indexing all

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
        Map<String, Object> toIndex = new HashMap<>();

        // Standard fields also included:
        toIndex.put("contig", variant.getContig());
        toIndex.put("start", variant.getStart());
        toIndex.put("end", variant.getEnd());
        // TODO: maybe others?


        for (String infoField : infoFieldsToIndex) {
            if (variant.hasAttribute(infoField) && variant.getAttribute(infoField) != null) {
                toIndex.put(infoField, variant.getAttribute(infoField));
            } else {
                throw new GATKException("Non-existent INFO field key: " + infoField);
            }
        }

        for (InfoFieldAnnotation a : annotationsToUse) {
            Map<String, Object> infoField = a.annotate(referenceContext, variant, null);

            for (var entry : infoField.entrySet()) {
                if(toIndex.containsKey(entry.getKey())) {
                    throw new GATKException("Duplicate INFO field key: " + entry.getKey());
                }
            }

            toIndex.putAll(infoField);
        }

        sites++;

        Document doc = new Document();

        for (var entry : toIndex.entrySet()) {
            if(entry.getValue().getClass() == String.class) {
                doc.add(new TextField(entry.getKey(), (String) entry.getValue(), Field.Store.YES));
            } else if(entry.getValue().getClass() == Integer.class) {
                doc.add(new IntPoint(entry.getKey(), (int) entry.getValue()));
                doc.add(new StoredField(entry.getKey(), (int) entry.getValue()));
            } else if(entry.getValue().getClass() == Float.class) {
                doc.add(new FloatPoint(entry.getKey(), (float) entry.getValue()));
                doc.add(new StoredField(entry.getKey(), (float) entry.getValue()));
            }
        }

        try {
            writer.addDocument(doc);
        } catch(IOException e) {
            throw new GATKException(e.getMessage(), e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Indexing complete, total sites indexed: " + sites);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        try {
            writer.close();
        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }
    }
}
