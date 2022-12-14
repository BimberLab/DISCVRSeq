package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.github.discvrseq.walkers.annotator.DiscvrVariantAnnotator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.*;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.store.FSDirectory;
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

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;

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
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", optional = true)
    File outDir;

    @Argument(doc="Info fields to index", fullName = "info-field", shortName = "IF", optional = true)
    List<String> infoFieldsToIndex;

    @Argument(doc="Annotations to apply", fullName = "annotation-name", shortName = "AN", optional = true)
    List<String> annotationClassNames;

    @Argument(fullName = "threads", doc="The number of threads to use.", optional=true)
    public int threads = 1;

    private IndexWriter writer = null;

    private FSDirectory index = null;

    private StandardAnalyzer analyzer = null;

    private VCFHeader header;

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
        analyzer = new StandardAnalyzer();

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

        header = (VCFHeader) getHeaderForFeatures(getDrivingVariantsFeatureInput());

        for(String field : infoFieldsToIndex) {
           if(!header.hasInfoLine(field)) {
                throw new GATKException("Non-existent INFO field key: " + field);
           }
        }

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

        dict = getBestAvailableSequenceDictionary();

        if (threads > 1) {
            executor = Executors.newScheduledThreadPool(threads);
        }
    }

    private ScheduledExecutorService executor = null;

    private SAMSequenceDictionary dict = null;

    private long sites = 0;

    private int getGenomicPosition(String contig, int start){
        int pos = start;

        // Add the length of each prior contig:
        for (SAMSequenceRecord sr : dict.getSequences()) {
            if (sr.getSequenceName().equals(contig)) {
                break;
            }

            pos += sr.getSequenceLength();
        }

        return pos;
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (executor != null) {
            try {
                executor.invokeAll(Collections.singletonList(new ApplyRunner(variant, referenceContext)));
            } catch (InterruptedException e) {
                throw new IllegalStateException("Error running VariantQC", e);
            }
        } else {
            new ApplyRunner(variant, referenceContext).call();
        }

        sites++;
    }

    private class ApplyRunner implements Callable<Boolean> {
        final VariantContext variant;
        final ReferenceContext referenceContext;

        public ApplyRunner(VariantContext variant, ReferenceContext referenceContext) {
            this.variant = variant;
            this.referenceContext = referenceContext;
        }

        @Override
        public Boolean call() {
            Map<String, List<Object>> toIndex = new HashMap<>();

            for (String infoField : infoFieldsToIndex) {
                if (variant.hasAttribute(infoField) && variant.getAttribute(infoField) != null) {
                    toIndex.put(infoField, variant.getAttributeAsList(infoField));
                }
            }

            for (InfoFieldAnnotation a : annotationsToUse) {
                Map<String, Object> infoField = a.annotate(referenceContext, variant, null);

                for (Map.Entry<String, Object> entry : infoField.entrySet()) {
                    if (toIndex.containsKey(entry.getKey())) {
                        throw new GATKException("Duplicate INFO field key: " + entry.getKey());
                    }

                    toIndex.put(entry.getKey(), Collections.singletonList(entry.getValue()));
                }
            }

            Document doc = new Document();

            // Add standard fields
            doc.add(new TextField("contig", variant.getContig(), Field.Store.YES));

            doc.add(new IntPoint("start", variant.getStart()));
            doc.add(new StoredField("start", variant.getStart()));

            doc.add(new IntPoint("end", variant.getEnd()));
            doc.add(new StoredField("end", variant.getEnd()));

            doc.add(new StoredField("genomicPosition", getGenomicPosition(variant.getContig(), variant.getStart())));

            addFieldsToDocument(doc, header, toIndex);

            try {
                writer.addDocument(doc);
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            return true;
        }
    }

    public static void addFieldsToDocument(Document doc, VCFHeader header, Map<String, List<Object>> toIndex) {
        for (Map.Entry<String, List<Object>> entry : toIndex.entrySet()) {
            for (Object value : entry.getValue()) {
                VCFInfoHeaderLine info = header.getInfoHeaderLine(entry.getKey());
                VCFHeaderLineType datatype = info == null ? VCFHeaderLineType.String : info.getType();
                addFieldToDocument(doc, datatype, entry.getKey(), value);
            }
        }
    }

    public static void addFieldToDocument(Document doc, VCFHeaderLineType variantHeaderLineType, String key, Object value) {
        switch(variantHeaderLineType) {
            case Character:
                doc.add(new StringField(key, (String) value, Field.Store.YES));
                break;
            case Flag:
                doc.add(new StringField(key, (String) value, Field.Store.YES));
                break;
            case Float:
                doc.add(new FloatPoint(key, Float.parseFloat(value.toString())));
                doc.add(new StoredField(key, Float.parseFloat(value.toString())));
                break;
            case Integer:
                doc.add(new IntPoint(key, Integer.parseInt(value.toString())));
                doc.add(new StoredField(key, Integer.parseInt(value.toString())));
                break;
            case String:
                doc.add(new TextField(key, (String) value, Field.Store.YES));
                break;
            }
    }


    @Override
    public Object onTraversalSuccess() {
        logger.info("Indexing complete, total sites indexed: " + sites);

        if (executor != null) {
            executor.shutdown();
        }

        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        try {
            if (writer != null) {
                writer.close();
            }
        } catch (Throwable e) {
            logger.error("Unable to close writer", e);
        }

        try {
            if (index != null)  {
                index.close();
            }
        } catch (Throwable e) {
            logger.error("Unable to close index directory", e);
        }

        try {
            if (analyzer != null) {
                analyzer.close();
            }
        } catch (Throwable e) {
            logger.error("Unable to close analyzer", e);
        }
    }
}