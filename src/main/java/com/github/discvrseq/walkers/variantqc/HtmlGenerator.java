package com.github.discvrseq.walkers.variantqc;

import com.github.discvrseq.Main;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import javax.annotation.Nullable;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;

/**
 * Created by bimber on 5/18/2017.
 */
public class HtmlGenerator {
    public static final String[] JS_SCRIPTS = new String[]{
            "templates/assets/js/packages/jquery-3.7.0.min.js",
            "templates/assets/js/packages/jquery-ui.min.js",
            "templates/assets/js/packages/numeral.min.js",
            //"https://cdnjs.cloudflare.com/ajax/libs/numeral.js/2.0.6/numeral.min.js",
            "templates/assets/js/packages/bootstrap.min.js",
            "templates/assets/js/packages/highcharts.js",
            "templates/assets/js/packages/highcharts.exporting.js",
            "templates/assets/js/packages/highcharts.heatmap.js",
            "templates/assets/js/packages/highcharts.offline-exporting.js",
            "templates/assets/js/packages/jquery.tablesorter.min.js",
            "templates/assets/js/packages/clipboard.min.js",
            "templates/assets/js/packages/FileSaver.min.js",
            //"https://cdnjs.cloudflare.com/ajax/libs/chroma-js/1.3.3/chroma.min.js",
            "templates/assets/js/packages/chroma.min.js",
            "templates/assets/js/packages/lz-string.min.js",
            "templates/assets/js/summaryTable.js",
            "templates/assets/js/shim.js"
    };

    public static final String[] JS_SCRIPTS2 = new String[]{
            "templates/assets/js/multiqc_tables.js",
            "templates/assets/js/multiqc_toolbox.js",
            "templates/assets/js/multiqc.js",
            "templates/assets/js/multiqc_plotting.js"
    };

    public static final String[] CSS_FILES = new String[]{
            "templates/assets/css/bootstrap.min.css",
            "templates/assets/css/default_multiqc.css",
            "templates/assets/css/font.css"
    };

    public HtmlGenerator() {

    }

    //separated for testing purposes
    public static void printStaticContent(PrintWriter out) throws IOException {
        //append header
        Resource header = new Resource("templates/template1.html", VariantQC.class);
        IOUtils.copy(header.getResourceContentsAsStream(), out, StandardCharsets.UTF_8);

        //scripts:
        for (String script : CSS_FILES){
            Resource r = new Resource(script, VariantQC.class);
            out.println("<style>");
            IOUtils.copy(r.getResourceContentsAsStream(), out, StandardCharsets.UTF_8);
            out.println("</style>");
        }

        for (String script : JS_SCRIPTS){
            appendScript(script, out);
        }
    }

    private static final String DEV = "*DevelopmentVersion*";

    private String getVersion() {
        return Main.class.getPackage().getImplementationVersion() == null ? DEV : Main.class.getPackage().getImplementationVersion();
    }

    private boolean isDevOrTest() {
        return DEV.equals(getVersion());
    }

    public void generateHtml(Collection<SectionJsonDescriptor> translatorList, PrintWriter out, @Nullable PrintWriter jsonWriter) throws IOException {

        printStaticContent(out);

        //config:
        out.println("<script type=\"text/javascript\">");
        out.println("mqc_plots = {};");
        out.println("num_datasets_plot_limit = 50;");
        out.println("$(function() {");

        //NOTE: if this is a dev or test build, use a static value (instead of date), so we can compare against a static output file
        String dateStr = new SimpleDateFormat("yyyy-MM-dd HH:mm").format(new Date());
        out.println("$('#dateTime').html('Generated on: " + (isDevOrTest() ? "*Timestamp*" : dateStr) + "');");
        out.println("processPlots({");
        out.println("sections:");

        JsonArray arr = new JsonArray();
        for (SectionJsonDescriptor t : translatorList){
            arr.add(t.getConfig());
        }
        out.println(arr.toString());
        out.println("});");
        out.println("});");
        out.println("</script>");

        for (String script : JS_SCRIPTS2){
            appendScript(script, out);
        }

        //append header
        Resource header2 = new Resource("templates/template2.html", VariantQC.class);

        try (StringWriter writer = new StringWriter()) {
            IOUtils.copy(header2.getResourceContentsAsStream(), writer, StandardCharsets.UTF_8);

            String toWrite = writer.getBuffer().toString().replaceAll("\\{version}", getVersion());
            IOUtils.write(toWrite, out);
        }

        //write raw data
        if (jsonWriter != null) {
            JsonObject json = new JsonObject();
            arr.forEach(e -> {
                JsonObject j = e.getAsJsonObject();
                String section = j.get("label").getAsString();
                JsonObject sectionJson = json.has(section) ? json.get(section).getAsJsonObject() : new JsonObject();

                JsonArray reports = j.get("reports").getAsJsonArray();
                reports.forEach( report -> {
                    JsonObject r = report.getAsJsonObject();
                    String reportName = r.get("label").getAsString();
                    if (sectionJson.has(reportName)){
                        throw new RuntimeException("Duplicate report name: " + reportName);
                    }
                    sectionJson.add(reportName, r);
                });

                json.add(section, sectionJson);
            });

            Gson gson = new GsonBuilder().setPrettyPrinting().create();
            jsonWriter.write(gson.toJson(json));
        }

    }

    private static void appendScript(String script, PrintWriter out) throws IOException{
        Resource r = new Resource(script, VariantQC.class);
        out.println("<script type=\"text/javascript\">");
        IOUtils.copy(r.getResourceContentsAsStream(), out, StandardCharsets.UTF_8);
        out.println("</script>");
    }
}
