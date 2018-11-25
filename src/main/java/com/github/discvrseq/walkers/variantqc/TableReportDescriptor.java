package com.github.discvrseq.walkers.variantqc;

import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.hellbender.utils.report.GATKReportColumn;
import org.broadinstitute.hellbender.utils.report.GATKReportDataType;

import java.util.*;


/**
 * Created by bimber on 5/22/2017.
 */
public class TableReportDescriptor extends ReportDescriptor {
    private Set<String> skippedColNames = new HashSet<>();

    public TableReportDescriptor(String label, String evaluatorModuleName, Collection<String> skippedSamples) {
        super(label, SectionJsonDescriptor.PlotType.data_table, evaluatorModuleName);
        skippedColNames.add(evaluatorModuleName);
        skippedColNames.add("EvalFeatureInput");
        skippedColNames.add("CompFeatureInput");
        if (skippedSamples != null) {
            this.skippedSamples.addAll(skippedSamples);
        }
    }

    public TableReportDescriptor(String label, String evaluatorModuleName) {
        this(label, evaluatorModuleName, null);
    }

    public static TableReportDescriptor getCountVariantsTable(boolean skipAll) {
        TableReportDescriptor ret = new TableReportDescriptor("Variant Summary", "CountVariants", skipAll ? Arrays.asList("all") : null);

        //JsonObject myColJson = new JsonObject();
        //myColJson.addProperty("dmin", 0);
        //myColJson.addProperty("dmax", 1.0);
        //ret.addColumnInfo("myColumn", myColJson);

        return ret;
    }

    public static TableReportDescriptor getIndelTable() {
        TableReportDescriptor ret = new TableReportDescriptor("Indel Summary", "IndelSummary", Arrays.asList("all"));
        ret.skippedColNames.add("n_indels_matching_gold_standard");
        ret.skippedColNames.add("gold_standard_matching_rate");

        return ret;
    }

    @Override
    public JsonObject getReportJson(String sectionTitle) {
        JsonObject ret = new JsonObject();
        ret.addProperty("label", label);

        JsonObject dataObj = new JsonObject();
        ret.add("data", dataObj);

        dataObj.addProperty("plot_type", plotType.name());

        dataObj.add("samples", getSampleNames());//Ordering of sample names must correspond with dataset order

        JsonArray datasetsJson = new JsonArray();
        int rowIdx = 0;
        while (rowIdx < table.getNumRows()) {
            List<Object> rowList = new ArrayList<>();
            String sampleName = getSampleNameForRow(rowIdx);
            if (skippedSamples.contains(sampleName)) {
                continue;
            }

            for (GATKReportColumn col : table.getColumnInfo()) {
                if (skippedColNames.contains(col.getColumnName())) {
                    continue;
                }

                rowList.add(table.get(rowIdx, col.getColumnName()));
            }

            datasetsJson.add(new GsonBuilder().create().toJsonTree(rowList).getAsJsonArray());
        }
        dataObj.add("datasets", datasetsJson);

        dataObj.add("columns", new JsonArray());
        for (GATKReportColumn col : table.getColumnInfo()) {
            if (skippedColNames.contains(col.getColumnName())) {
                continue;
            }

            JsonObject colJson = new JsonObject();
            colJson.addProperty("name", col.getColumnName());
            colJson.addProperty("label", descriptionMap.containsKey(col.getColumnName()) ? descriptionMap.get(col.getColumnName()) : col.getColumnName());

            if (col.getDataType() == GATKReportDataType.Decimal) {
                //TODO: look into format strings supporting more than 6 decimals
                //colJson.addProperty("formatString", "0.0[0000]");

                flagValueByTwoStandardDeviation(colJson, col.getColumnName());
                inferMinMax(colJson, col.getColumnName());

            } else if (col.getDataType() == GATKReportDataType.Integer) {
                colJson.addProperty("formatString", "0,0");
                flagValueByTwoStandardDeviation(colJson, col.getColumnName());
                inferMinMax(colJson, col.getColumnName());
            }

            //allow upstream code to supply custom config
            if (columnInfoMap.containsKey(col.getColumnName())) {
                for (Map.Entry<String, JsonElement> e : columnInfoMap.get(col.getColumnName()).entrySet()) {
                    colJson.add(e.getKey(), e.getValue());
                }
            }

            dataObj.getAsJsonArray("columns").add(colJson);
        }

        return ret;

    }

    private void inferMinMax(JsonObject colJson, String colName) {
        List<Double> rowValuesList = new ArrayList<>();
        int rowIdx = 0;
        while (rowIdx < table.getNumRows()) {
            if (skippedSamples.contains(getSampleNameForRow(rowIdx))) {
                continue;
            }

            rowValuesList.add(NumberUtils.createNumber(table.get(rowIdx, colName).toString()).doubleValue());
        }

        Double min = Collections.min(rowValuesList) - Collections.min(rowValuesList) * 0.1;
        Double max = Collections.max(rowValuesList) + Collections.max(rowValuesList) * 0.1;
        colJson.addProperty("dmin", min);
        colJson.addProperty("dmax", max == 0 ? 1 : max);
    }

    private void flagValueByTwoStandardDeviation(JsonObject colJson, String colName){
        DescriptiveStatistics stats = new DescriptiveStatistics();
        int rowIdx = 0;
        while (rowIdx < table.getNumRows()) {
            if (skippedSamples.contains(getSampleNameForRow(rowIdx))) {
                continue;
            }

            stats.addValue(NumberUtils.createBigDecimal(table.get(rowIdx, colName).toString()).doubleValue());
        }

        Double sd = stats.getStandardDeviation();
        Double mean = stats.getMean();
        Double aboveTwoSd = mean + (2.0 * sd);
        Double belowTwoSd = mean - (2.0 * sd);
        colJson.addProperty("flagAbove", aboveTwoSd);
        colJson.addProperty("flagBelow", belowTwoSd);
    }
}