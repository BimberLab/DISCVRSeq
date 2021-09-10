package com.github.discvrseq.walkers.variantqc;

import com.google.gson.JsonObject;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.utils.report.GATKReportColumn;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.samples.SampleDB;

import javax.annotation.Nullable;
import java.util.*;

/**
 * Created by bimber on 5/22/2017.
 */
abstract class ReportDescriptor {
    protected final String sectionLabel;
    protected final String reportLabel;
    protected final SectionJsonDescriptor.PlotType plotType;
    protected final String evaluatorModuleName;
    protected Map<String, JsonObject> columnInfoMap;
    protected SectionJsonDescriptor sectionConfig;
    protected GATKReportTable table;
    protected Set<String> skippedSamples = new HashSet<>();
    protected Map<String, String> descriptionMap;
    protected GATKReportTableTransformer transformer;
    protected boolean isMultiVcf;
    protected Comparator<String> comparator = Comparator.naturalOrder(); //new NaturalSortComparator();

    protected ReportDescriptor(String reportLabel, String sectionLabel, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName, @Nullable GATKReportTableTransformer transformer, boolean isMultiVcf) {
        this.sectionLabel = sectionLabel;
        this.reportLabel = reportLabel;
        this.plotType = plotType;
        this.evaluatorModuleName = evaluatorModuleName;
        this.columnInfoMap = new HashMap<>();
        this.transformer = transformer;
        this.isMultiVcf = isMultiVcf;
    }

    protected void setComparator(Comparator<String> comparator) {
        this.comparator = comparator;
    }

    public String getEvaluatorModuleName() {
        return evaluatorModuleName;
    }

    public void addColumnInfo(String colName, JsonObject columnInfo) {
        this.columnInfoMap.put(colName, columnInfo);
    }

    public void bindSection(SectionJsonDescriptor sectionConfig, GATKReportTable table, Map<String, String> descriptionMap, SampleDB sampleDB){
        if (transformer != null){
            table = transformer.transform(table, sampleDB);
        }

        this.sectionConfig = sectionConfig;
        this.table = table;
        this.descriptionMap = descriptionMap;
    }

    abstract JsonObject getReportJson(String sectionTitle);

    protected List<String> columnsInSampleName = null;

    /**
     * This provides the opportunity for subclasses to supply custom logic to parse the sampleName from rows.
     * This can return null, in which case that row will be ignored, for example if specific states should not be included.
     * @param rowIdx The row index to parse
     * @return the parsed sample name
     */
    protected String getSampleNameForRow(int rowIdx){
        if (columnsInSampleName == null){
            List<String> ret = new ArrayList<>();
            for (GATKReportColumn col : table.getColumnInfo()){
                if (sectionConfig.stratifications.contains(col.getColumnName()) || (isMultiVcf && "EvalFeatureInput".equals(col.getColumnName()))){
                    ret.add(col.getColumnName());
                }
            }

            // Ensure this column is first:
            if (isMultiVcf && ret.contains("EvalFeatureInput")) {
                ret.remove("EvalFeatureInput");
                ret.add(0, "EvalFeatureInput");
            }
            columnsInSampleName = ret;
        }

        List<String> tokens = new ArrayList<>();
        for (String colName : columnsInSampleName){
            tokens.add(table.get(rowIdx, colName).toString());
        }

        return StringUtil.join(" / ", tokens);
    }

    protected boolean shouldSkipRow(int rowIdx) {
        if (skippedSamples.isEmpty()) {
            return false;
        }

        String sn = getSampleNameForRow(rowIdx);
        if (sn != null && skippedSamples.contains(sn)) {
            return true;
        }

        //Note: check both the Sample Name column (which will be a concatenation of sample and stratifiers, and also the simple 'Sample' column, if present.
        if (hasSampleColumn()) {
            sn = (String) table.get(rowIdx, "Sample");
            if (sn != null && skippedSamples.contains(sn)) {
                return true;
            }
        }

        return false;
    }

    protected Boolean hasSampleColumn = null;

    private boolean hasSampleColumn() {
        if (hasSampleColumn == null) {
            hasSampleColumn = false;
            for (GATKReportColumn col : table.getColumnInfo()){
                if ("Sample".equals(col.getColumnName())){
                    hasSampleColumn = true;
                    break;
                }
            }
        }

        return hasSampleColumn;
    }

    protected Map<String, Integer> getSampleToRowIdx(){
        Map<String, Integer> sampleToRowIdx = new TreeMap<>(comparator);
        for (int rowIdx = 0;rowIdx<table.getNumRows();rowIdx++){
            if (!shouldSkipRow(rowIdx)){
                sampleToRowIdx.put(getSampleNameForRow(rowIdx), rowIdx);
            }
        }

        return sampleToRowIdx;
    }
}
