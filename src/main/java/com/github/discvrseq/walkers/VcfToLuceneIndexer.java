package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.ICSVWriter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.*;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.search.Sort;
import org.apache.lucene.search.SortField;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.BytesRef;
import org.apache.lucene.util.NumericUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;

import javax.annotation.Nullable;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

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
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool accepts a VCF, iterates the variants, and writes the results to a lucene search index",
        oneLineSummary = "Index a VCF using Apache Lucene",
        programGroup = DiscvrSeqProgramGroup.class
)
public class VcfToLuceneIndexer extends VariantWalker {
    public static final int MAX_VALUES_TO_PRINT = 200;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", optional = true)
    GATKPath outDir;

    @Argument(fullName = "index-stats", doc = "A file where a TSV of summary information about each indexed field will be written, including min/max value and a list of unique values", optional = true)
    GATKPath indexStatsPath;

    @Argument(doc="Info fields to index", fullName = "info-field", shortName = "IF", optional = true)
    List<String> infoFieldsToIndex;

    @Argument(fullName = "threads", doc="The number of threads to use.", optional=true)
    public int threads = 1;

    @Argument(fullName = "allow-missing-fields", doc="If true, the tool will warn, rather than fail, if a non-existent --info-field is requested.", optional=true)
    public boolean allowMissingFields = false;

    @Argument(fullName= "validation-stringency", doc = "The level of validation, either LENIENT or STRICT", common = true, optional = true)
    protected ValidationStringency stringency = ValidationStringency.STRICT;

    private IndexWriter writer = null;

    private final IndexStats stats = new IndexStats();

    private FSDirectory index = null;

    private StandardAnalyzer analyzer = null;

    @Override
    public boolean useVariantAnnotations() {
        return true;
    }

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
        config.setIndexSort(new Sort(new SortField("genomicPosition_sort", SortField.Type.LONG, false)));

        try {
            writer = new IndexWriter(index, config);
        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        VCFHeader header = (VCFHeader) getHeaderForFeatures(getDrivingVariantsFeatureInput());

        List<String> missing = new ArrayList<>();
        for (String field : infoFieldsToIndex) {
            if (!header.hasInfoLine(field)) {
                if (allowMissingFields) {
                    missing.add(field);
                }
                else {
                    throw new GATKException("Non-existent INFO field key: " + field);
                }
            }
            else {
                stats.addField(header.getInfoHeaderLine(field));
            }
        }

        if (!missing.isEmpty()) {
            logger.warn("The following fields were requested but not present: " + StringUtils.join(missing, ","));
        }

        if (threads > 1) {
            executor = Executors.newScheduledThreadPool(threads);
        }

        // populate map. Add the length of each prior contig:
        SAMSequenceDictionary dict = getBestAvailableSequenceDictionary();
        for (SAMSequenceRecord sr1 : dict.getSequences()) {
            long offset = 0;
            for (SAMSequenceRecord sr : dict.getSequences()) {
                if (sr.getSequenceName().equals(sr1.getSequenceName())) {
                    break;
                }

                offset += sr.getSequenceLength();
            }
            genomicPositionMap.put(sr1.getSequenceName(), offset);
        }
    }

    private ScheduledExecutorService executor = null;

    private long sites = 0;

    private final Map<String, Long> genomicPositionMap = new HashMap<>();

    private long getGenomicPosition(String contig, int start){
        return start + genomicPositionMap.get(contig);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (executor != null) {
            try {
                executor.invokeAll(Collections.singletonList(new ApplyRunner(variant, referenceContext)));
            } catch (InterruptedException e) {
                throw new IllegalStateException("Error running VcfToLuceneIndexer", e);
            }
        } else {
            new ApplyRunner(variant, referenceContext).call();
        }

        sites++;
    }

    public class ApplyRunner implements Callable<Boolean> {
        final VariantContext variant;
        final ReferenceContext referenceContext;

        public ApplyRunner(VariantContext variant, ReferenceContext referenceContext) {
            this.variant = variant;
            this.referenceContext = referenceContext;
        }

        @Override
        public Boolean call() {
            final VCFHeader header = getHeaderForVariants();

            // Index each ALT by itself:
            for (Allele alt : variant.getAlternateAlleles()) {
                Document doc = new Document();
                int alleleIdx = variant.getAlleleIndex(alt); // includes REF
                int altAlleleIndex = alleleIdx - 1;

                for (String infoField : infoFieldsToIndex) {
                    if (variant.hasAttribute(infoField) && variant.getAttribute(infoField) != null) {
                        final VCFInfoHeaderLine line = header.getInfoHeaderLine(infoField);
                        final VCFHeaderLineType datatype = line.getType();

                        if (line.getCountType() == VCFHeaderLineCount.A) {
                            List<Object> vals = variant.getAttributeAsList(infoField);
                            if (vals.size() != variant.getAlternateAlleles().size()) {
                                String msg = "Incorrect number of annotations for " + infoField + ". Expected: " + variant.getAlternateAlleles().size() + ", found: " + vals.size() + ". Value: " + variant.getAttribute(infoField) + ", at " + variant.toStringWithoutGenotypes();
                                if (stringency == ValidationStringency.STRICT) {
                                    throw new GATKException(msg);
                                }
                                else {
                                    logger.warn(msg + ", skipping");
                                    continue;
                                }
                            }

                            Object val = vals.get(altAlleleIndex);
                            if (val != null) {
                                addFieldToDocument(doc, datatype, infoField, val);
                            }
                        }
                        else if (line.getCountType() == VCFHeaderLineCount.R) {
                            List<Object> vals = variant.getAttributeAsList(infoField);
                            if (vals.size() != variant.getNAlleles()) {
                                String msg = "Incorrect number of annotations for " + infoField + ". Expected: " + variant.getNAlleles() + ", found: " + vals.size() + ". Value: " + variant.getAttribute(infoField) + ", at " + variant.toStringWithoutGenotypes();
                                if (stringency == ValidationStringency.STRICT) {
                                    throw new GATKException(msg);
                                }
                                else {
                                    logger.warn(msg + ", skipping");
                                    continue;
                                }
                            }

                            Object val = vals.get(alleleIdx);
                            if (val != null) {
                                addFieldToDocument(doc, datatype, infoField, val);
                            }
                        }
                        else if (line.getCountType() == VCFHeaderLineCount.INTEGER || line.getCountType() == VCFHeaderLineCount.UNBOUNDED) {
                            Object val = variant.getAttribute(infoField);
                            addFieldToDocument(doc, datatype, infoField, val);
                        }
                    }
                }

                // Add standard fields
                doc.add(new TextField("contig", variant.getContig(), Field.Store.YES));
                doc.add(new SortedDocValuesField("contig_sort", new BytesRef(variant.getContig())));

                doc.add(new TextField("ref", variant.getReference().getDisplayString(), Field.Store.YES));
                doc.add(new SortedDocValuesField("ref_sort", new BytesRef(variant.getReference().getDisplayString())));

                doc.add(new TextField("alt", alt.getDisplayString(), Field.Store.YES));
                doc.add(new SortedDocValuesField("alt_sort", new BytesRef(alt.getDisplayString())));

                final long genomicPositionStart = getGenomicPosition(variant.getContig(), variant.getStart());
                doc.add(new IntPoint("start", variant.getStart()));
                doc.add(new StoredField("start", variant.getStart()));
                doc.add(new NumericDocValuesField("start_sort", genomicPositionStart));

                final long genomicPositionEnd = getGenomicPosition(variant.getContig(), variant.getEnd());
                doc.add(new IntPoint("end", variant.getEnd()));
                doc.add(new StoredField("end", variant.getEnd()));
                doc.add(new NumericDocValuesField("end_sort", genomicPositionEnd));

                doc.add(new LongPoint("genomicPosition", genomicPositionStart));
                doc.add(new StoredField("genomicPosition", genomicPositionStart));
                doc.add(new NumericDocValuesField("genomicPosition_sort", genomicPositionStart));

                if (variant.hasGenotypes()) {
                    AtomicReference<String> docValue = new AtomicReference<>(null);

                    variant.getGenotypes().stream().filter(g -> !g.isFiltered() && !g.isNoCall() && g.getAlleles().contains(alt)).map(Genotype::getSampleName).sorted().forEach(sample -> {
                        doc.add(new TextField("variableSamples", sample, Field.Store.YES));

                        if (docValue.get() == null) {
                            docValue.set(sample);
                        }
                    });

                    if (docValue.get() != null) {
                        doc.add(new SortedDocValuesField("variableSamples_sort", new BytesRef(docValue.get())));
                        docValue.set(null);
                    }

                    variant.getGenotypes().stream().filter(g -> !g.isFiltered() && !g.isNoCall() && g.getAlleles().contains(alt) && g.isHomVar()).map(Genotype::getSampleName).sorted().forEach(sample -> {
                        doc.add(new TextField("homozygousVarSamples", sample, Field.Store.YES));

                        if (docValue.get() == null) {
                            docValue.set(sample);
                        }
                    });

                    if (docValue.get() != null) {
                        doc.add(new SortedDocValuesField("homozygousVarSamples_sort", new BytesRef(docValue.get())));
                        docValue.set(null);
                    }

                    long nHet = variant.getGenotypes().stream().filter(g -> !g.isFiltered() && !g.isNoCall() && g.getAlleles().contains(alt) && g.isHet()).count();
                    doc.add(new IntPoint("nHet", (int)nHet));
                    doc.add(new StoredField("nHet", (int)nHet));
                    doc.add(new NumericDocValuesField("nHet_sort", (int)nHet));

                    long nHomVar = variant.getGenotypes().stream().filter(g -> !g.isFiltered() && !g.isNoCall() && g.getAlleles().contains(alt) && g.isHomVar()).count();
                    doc.add(new IntPoint("nHomVar", (int)nHomVar));
                    doc.add(new StoredField("nHomVar", (int)nHomVar));
                    doc.add(new NumericDocValuesField("nHomVar_sort", (int)nHomVar));

                    long nCalled = variant.getGenotypes().stream().filter(g -> !g.isFiltered() && !g.isNoCall()).count();
                    doc.add(new IntPoint("nCalled", (int)nCalled));
                    doc.add(new StoredField("nCalled", (int)nCalled));
                    doc.add(new NumericDocValuesField("nCalled_sort", (int)nCalled));

                    float fractionHet = (float) nHet / (float) (nHet + nHomVar);
                    doc.add(new DoublePoint("fractionHet", fractionHet));
                    doc.add(new StoredField("fractionHet", fractionHet));
                    doc.add(new NumericDocValuesField("fractionHet_sort", NumericUtils.doubleToSortableLong(fractionHet)));
                }

                try {
                    writer.addDocument(doc);
                } catch (IOException e) {
                    throw new GATKException(e.getMessage(), e);
                }
            }

            return true;
        }
    }

    private <T> @Nullable Collection<T> attemptToFixNumericValue(String key, Object value, Class<T> clazz) {
        // NOTE: there are situations where a numeric value can have duplicate values for a given variant/allele
        // This is sort of a hack, but in this situation we will just double-index them:
        String valueStr = value.toString();
        if (!NumberUtils.isCreatable(valueStr)) {
            if (valueStr.contains("|")) {
                try {
                    if (clazz == Double.class) {
                        return Arrays.stream(valueStr.split("\\|")).filter(x -> !x.isEmpty()).map(Double::parseDouble).map(clazz::cast).collect(Collectors.toSet());
                    } else if (clazz == Integer.class) {
                        return Arrays.stream(valueStr.split("\\|")).filter(x -> !x.isEmpty()).map(Integer::parseInt).map(clazz::cast).collect(Collectors.toSet());
                    }

                    // should never reach this
                    throw new GATKException("Attempted to parse a numeric value into something other than Double/Integer. This should never occur.");
                } catch (Exception e) {
                    possiblyReportBadValue(e, key, valueStr);
                    return null;
                }
            }
            else {
                possiblyReportBadValue(null, key, valueStr);
                return null;
            }
        }
        else {
            try {
                if (clazz == Double.class) {
                    return Collections.singleton(clazz.cast(Double.parseDouble(valueStr)));
                }
                else if (clazz == Integer.class) {
                    return Collections.singleton(clazz.cast(Integer.parseInt(valueStr)));
                }

                throw new GATKException("Attempted to parse a numeric value into something other than Double/Integer. This should never occur.");
            }
            catch (Exception e) {
                possiblyReportBadValue(e, key, valueStr);
                return null;
            }
        }
    }

    private final Set<String> keysWithErrors = new HashSet<>();

    private void possiblyReportBadValue(@Nullable Exception e, String key, Object fieldValue) {
        String message = "Unable to parse field: " + key + ", was: <" + fieldValue + ">. " + (e == null ? "" : e.getMessage());
        if (stringency == ValidationStringency.STRICT) {
            throw new GATKException(message);
        }
        else {
            if (!keysWithErrors.contains(key)) {
                keysWithErrors.add(key);
                logger.warn(message);
            }
        }
    }

    synchronized private void addFieldToDocument(Document doc, VCFHeaderLineType variantHeaderLineType, String key, Object fieldValue) {
        try {
            stats.inspectValue(key, fieldValue);
        }
        catch (Exception e) {
            possiblyReportBadValue(e, key, fieldValue);
        }

        Collection<?> values = fieldValue instanceof Collection ? (Collection<?>) fieldValue : Collections.singleton(fieldValue);

        AtomicBoolean indexDocValue = new AtomicBoolean(true);

        values.forEach(value -> {
            if (value == null || "".equals(value) || VCFConstants.EMPTY_INFO_FIELD.equals(value)) {
                return;
            }

            try {
                switch (variantHeaderLineType) {
                    case Character -> {
                        doc.add(new StringField(key, String.valueOf(value), Field.Store.YES));

                        if (indexDocValue.get()) {
                            doc.add(new SortedDocValuesField(key + "_sort", new BytesRef(String.valueOf(value))));
                            indexDocValue.set(false);
                        }
                    }
                    case Flag -> {
                        int x = Boolean.parseBoolean(value.toString()) ? 1 : 0;
                        doc.add(new IntPoint(key, x));

                        if (indexDocValue.get()) {
                            doc.add(new NumericDocValuesField(key + "_sort", x));
                            indexDocValue.set(false);
                        }
                    }
                    case Float -> {
                        Collection<Double> parsedVals = attemptToFixNumericValue(key, value, Double.class);
                        if (parsedVals != null) {
                            AtomicReference<Double> docValue = new AtomicReference<>(null);

                            parsedVals.forEach(x -> {
                                doc.add(new DoublePoint(key, x));
                                doc.add(new StoredField(key, x));

                                if (docValue.get() == null) {
                                    docValue.set(x);
                                }
                            });

                            if (docValue.get() != null && indexDocValue.get()) {
                                doc.add(new NumericDocValuesField(key + "_sort", NumericUtils.doubleToSortableLong(docValue.get())));
                                indexDocValue.set(false);
                            }
                        }
                    }
                    case Integer -> {
                        Collection<Integer> parsedVals = attemptToFixNumericValue(key, value, Integer.class);
                        if (parsedVals != null) {
                            AtomicReference<Integer> docValue = new AtomicReference<>(null);
                            parsedVals.forEach(x -> {
                                doc.add(new IntPoint(key, x));
                                doc.add(new StoredField(key, x));

                                if (docValue.get() == null) {
                                    docValue.set(x);
                                }
                            });

                            if (docValue.get() != null && indexDocValue.get()) {
                                doc.add(new NumericDocValuesField(key + "_sort", docValue.get()));
                                indexDocValue.set(false);
                            }
                        }
                    }
                    case String -> {
                        doc.add(new TextField(key, String.valueOf(value), Field.Store.YES));

                        if (indexDocValue.get()) {
                            doc.add(new SortedDocValuesField(key +"_sort", new BytesRef(String.valueOf(value))));
                            indexDocValue.set(false);
                        }
                    }
                    default -> possiblyReportBadValue(new Exception("VCF header type was not expected: " + variantHeaderLineType.name()), key, value);
                }
            }
            catch (Exception e) {
                if (stringency == ValidationStringency.STRICT) {
                    throw e;
                }
                else {
                    logger.warn("Error parsing value: " + value + ", for key: " + key);
                }
            }
        });
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Indexing complete, total sites indexed: " + sites);

        if (executor != null) {
            executor.shutdown();
        }

        if (indexStatsPath != null) {
            try (ICSVWriter csvWriter = CsvUtils.getTsvWriter(indexStatsPath.toPath())) {
                csvWriter.writeNext(new String[]{"Key", "Type", "TotalIndexed", "ContainedMultiValuedRow", "MinVal", "MaxVal", "DistinctValues"});

                for (String key : stats.collectorMap.keySet()) {
                    IndexStats.Collector c = stats.collectorMap.get(key);

                    csvWriter.writeNext(c.getCsvRow(key));
                }
            }
            catch (IOException e) {
                throw new GATKException("Error generating index stats", e);
            }

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

    public static class IndexStats
    {
        private final Map<String, Collector> collectorMap = new HashMap<>();

        public void addField(VCFInfoHeaderLine line) {
            switch (line.getType()) {
                case Character, String, Flag -> collectorMap.put(line.getID(), new StringCollector());
                case Integer, Float -> collectorMap.put(line.getID(), new NumericCollector());
            }
        }

        public void inspectValue(String key, Object val) {
            collectorMap.get(key).inspect(val);
        }

        public abstract static class Collector {
            protected boolean containedMultiValue = false;
            protected long totalIndexed = 0L;

            public void inspect(Object val) {
                if (val == null || "".equals(val) || VCFConstants.EMPTY_INFO_FIELD.equals(val)) {
                    return;
                }

                if (val instanceof Collection<?> valCollection) {
                    containedMultiValue = true;
                    valCollection.forEach(this::inspectValue);
                }
                else {
                    inspectValue(val);
                }
            }

            abstract protected void inspectValue(Object val);

            abstract public String[] getCsvRow(String key);
        }

        private static class NumericCollector extends Collector {
            Double minVal = null;
            Double maxVal = null;

            @Override
            protected void inspectValue(Object val) {
                double d;
                if (val instanceof Number number) {
                    d = number.doubleValue();
                }
                else {
                    try {
                        d = Double.parseDouble(String.valueOf(val));
                    }
                    catch (Exception e) {
                        // Ignore, let reporting upstream handle this
                        return;
                    }
                }

                totalIndexed++;
                if (minVal == null || d < minVal) {
                    minVal = d;
                }

                if (maxVal == null || d > maxVal) {
                    maxVal = d;
                }
            }

            @Override
            public String[] getCsvRow(String key) {
                return new String[]{key, "Numeric", String.valueOf(totalIndexed), String.valueOf(containedMultiValue), minVal == null ? "" : String.valueOf(minVal), maxVal == null ? "" : String.valueOf(maxVal), ""};
            }
        }

        private static class StringCollector extends Collector {
            Set<Object> values = new HashSet<>();

            @Override
            protected void inspectValue(Object val) {
                if (val == null || "".equals(val) || VCFConstants.EMPTY_INFO_FIELD.equals(val)) {
                    return;
                }

                totalIndexed++;
                values.add(String.valueOf(val));
            }

            @Override
            public String[] getCsvRow(String key) {
                return new String[]{key, "String", String.valueOf(totalIndexed), String.valueOf(containedMultiValue), "", "", values.size() > MAX_VALUES_TO_PRINT ? "Too many: " + values.size() : StringUtils.join(values, ", ")};
            }
        }
    }
}