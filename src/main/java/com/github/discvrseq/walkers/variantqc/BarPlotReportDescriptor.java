package com.github.discvrseq.walkers.variantqc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;
import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.hellbender.utils.report.GATKReportColumn;
import org.broadinstitute.hellbender.utils.report.GATKReportDataType;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.util.*;


/**
 * Created by bimber on 5/22/2017.
 */
public class BarPlotReportDescriptor extends ReportDescriptor {
    private String[] columnsToPlot;
    private String yLabel;

    public BarPlotReportDescriptor(String plotTitle, String sectionLabel, boolean isMultiVcf, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName, String[] columnsToPlot, String yLabel, Collection<String> skippedSamples, PivotingTransformer transformer) {
        super(plotTitle, sectionLabel, plotType, evaluatorModuleName, transformer, isMultiVcf);
        this.columnsToPlot = columnsToPlot;
        this.yLabel = yLabel;
        if (skippedSamples != null){
            this.skippedSamples.addAll(skippedSamples);
        }
    }

    public List<String> getColumnsToPlot(GATKReportTable table){
        return Arrays.asList(columnsToPlot);
    }

    public static BarPlotReportDescriptor getVariantTypeBarPlot(String sectionLabel, boolean isMultiVcf) {
        return new BarPlotReportDescriptor("Variant Type", sectionLabel, isMultiVcf, SectionJsonDescriptor.PlotType.bar_graph, "CountVariants", new String[]{"nSNPs", "nMNPs", "nInsertions", "nDeletions", "nComplex", "nSymbolic", "nMixed"}, "# Variants", Arrays.asList("all"), null);
    }

    public static BarPlotReportDescriptor getSiteFilterTypeBarPlot(String sectionLabel, boolean isMultiVcf, PivotingTransformer transformer) {
        return new BarPlotReportDescriptor("Filter Type", sectionLabel, isMultiVcf, SectionJsonDescriptor.PlotType.bar_graph, "CountVariants", null, "# Variants", Arrays.asList("all"), transformer){
            @Override
            public List<String> getColumnsToPlot(GATKReportTable table){
                List<String> ret = new ArrayList<>();
                for (GATKReportColumn col : table.getColumnInfo()){
                    if (!sectionConfig.stratifications.contains(col.getColumnName()) && !evaluatorModuleName.equals(col.getColumnName())){
                        if (GATKReportDataType.Decimal == col.getDataType() || GATKReportDataType.Integer == col.getDataType())
                            ret.add(col.getColumnName());
                    }
                }

                return ret;
            }
        };
    }

    @Override
    public JsonObject getReportJson(String sectionTitle) {
        JsonObject ret = new JsonObject();
        ret.addProperty("label", reportLabel);

        JsonObject dataObj = new JsonObject();
        ret.add("data", dataObj);

        dataObj.addProperty("plot_type", plotType.name());

        dataObj.add("samples", new JsonArray());
        dataObj.getAsJsonArray("samples").add(getSampleNames());

        JsonArray datasetsJson = new JsonArray();

        for (String colName : getColumnsToPlot(table)) {
            JsonObject datasetJson = new JsonObject();
            datasetJson.addProperty("name", colName);

            JsonArray data = new JsonArray();
            for (int rowIdx = 0;rowIdx<table.getNumRows();rowIdx++){
                if (shouldSkipRow(rowIdx)){
                    continue;
                }

                data.add(new JsonPrimitive(NumberUtils.createNumber(table.get(rowIdx, colName).toString())));
            }
            datasetJson.add("data", data);

            datasetsJson.add(datasetJson);
        }
        dataObj.add("datasets", new JsonArray());
        dataObj.getAsJsonArray("datasets").add(datasetsJson);

        JsonObject configJson = new JsonObject();
        configJson.addProperty("ylab", this.yLabel);
        configJson.addProperty("title", reportLabel);

        dataObj.add("config", configJson);

        return ret;
    }
}
