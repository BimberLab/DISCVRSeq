package com.github.discvrseq.walkers.variantqc;

import com.github.discvrseq.util.NaturalSortComparator;
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
    private static final Comparator<String> NATURAL_SORT = new NaturalSortComparator();

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

    protected ReportDescriptor(String reportLabel, String sectionLabel, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName, @Nullable GATKReportTableTransformer transformer, boolean isMultiVcf) {
        this.sectionLabel = sectionLabel;
        this.reportLabel = reportLabel;
        this.plotType = plotType;
        this.evaluatorModuleName = evaluatorModuleName;
        this.columnInfoMap = new HashMap<>();
        this.transformer = transformer;
        this.isMultiVcf = isMultiVcf;
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

    protected List<String> getColumnsInSampleName()
    {
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

        return columnsInSampleName;
    }

    /**
     * This provides the opportunity for subclasses to supply custom logic to parse the sampleName from rows.
     * This can return null, in which case that row will be ignored, for example if specific states should not be included.
     * @param rowIdx The row index to parse
     * @return the parsed sample name
     */
    protected String getSampleNameForRow(int rowIdx){
        List<String> tokens = new ArrayList<>();
        for (String colName : getColumnsInSampleName()){
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

    protected Map<Integer, String> getSortedRows(){
        Map<Integer, List<String>> sortValues = new HashMap<>();
        for (int rowIdx = 0;rowIdx<table.getNumRows();rowIdx++){
            List<String> vals = new ArrayList<>();
            for (String col : getColumnsInSampleName()) {
                vals.add(String.valueOf(table.get(rowIdx, col)));
            }

            sortValues.put(rowIdx, vals);
        }

        List<Integer> rowIdxs = new ArrayList<>(sortValues.keySet());
        Collections.sort(rowIdxs, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                List<String> l1 = sortValues.get(o1);
                List<String> l2 = sortValues.get(o2);

                for (int idx=0;idx<l1.size();idx++) {
                    // special-case all to show first
                    if ("all".equalsIgnoreCase(l1.get(idx)) && !"all".equalsIgnoreCase(l2.get(idx))) {
                        return 1;
                    }

                    int sortVal = NATURAL_SORT.compare(l1.get(idx), l2.get(idx));
                    if (sortVal != 0) {
                        return sortVal;
                    }
                }

                return 0;
            }
        });

        Map<Integer, String> sampleToRowIdx = new LinkedHashMap<>();
        for (int rowIdx : rowIdxs){
            if (!shouldSkipRow(rowIdx)){
                sampleToRowIdx.put(rowIdx, getSampleNameForRow(rowIdx));
            }
        }

        return sampleToRowIdx;
    }
}
