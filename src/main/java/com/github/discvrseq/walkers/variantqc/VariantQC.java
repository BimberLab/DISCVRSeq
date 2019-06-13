package com.github.discvrseq.walkers.variantqc;

import au.com.bytecode.opencsv.CSVReader;
import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager.Stratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.AnalysisModuleScanner;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.report.GATKReportVersion;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.reflections.Reflections;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.*;
/**
 * VariantQC will generate a user-friendly HTML report that aggregates data from a VCF file by site, filter status, sample, and other stratifiers in order
 * to provide various summary tables and graphs.  In most cases the tables overlay bar graphs for numeric columns in order to help identify outliers.
 * This tool is analogous to FASTQC or MultiQC, which provide similar HTML summary reports for different types of sequence data.
 *
 * <h3>Usage examples:</h3>
 * <h4>Basic Usage, Without Pedigree:</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R Homo_sapiens_assembly38.fasta \
 *     -V input.vcf \
 *     -O output.html
 * </pre>
 *
 * <h4>Report With Pedigree (required for gender information to display):</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R Homo_sapiens_assembly38.fasta \
 *     -ped myPedigree.ped \
 *     -V input.vcf \
 *     -O output.html
 * </pre>

 * <h3>Example Report:</h3>
 * Below is an image of an example report, showing a VCF summarized by sample.  The left-hand navigation allows the user to toggle between different stratifications (entire VCF, sample, contig, etc.).  The report
 * contains a series of tables or graphs summarizing different aspects of the data.  A complete example HTML report can be <a href="resources/variantQCSampleReport.html">viewed here</a>
 * <p></p>
 * <img src="images/variantQCSampleReport1.png" style="width: 75%"/>
 *
 * <h3>Advanced Usage:</h3>
 *
 * TODO
 *
 * <h3>Authors:</h3>
 * <ul>
 *  <li>Ben Bimber, Ph.D.</li>
 *  <li>Melissa Yan</li>
 * </ul>
 *
 * <h3>Acknowledgements:</h3>
 * This tool internally uses <a href="https://github.com/broadinstitute/gatk">GATK4's VariantEval</a> to aggregate data, and borrows heavily from <a href="https://multiqc.info/">MultiQC</a> for the HTML in the final report.
 * We thank Philip Ewels for the use of MultiQC code, and Chris Norman of the Broad Institute for assistance with VariantEval.
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This will generate an HTML summary report for a given VCF, analogous to the report FastQC creates for raw read data.",
        oneLineSummary = "Generate HTML summary report for a VCF",
        programGroup = DiscvrSeqProgramGroup.class
)
public class VariantQC extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName = StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for identifying Mendelian violations", optional=true)
    private File pedigreeFile;

    @Argument(fullName = "pedigreeValidationType", shortName = "pedValidationType", doc="The strictness for validating the pedigree.  Can be either STRICT or SILENT.  Default is STRICT", optional=true)
    private PedigreeValidationType pedigreeValidationType = PedigreeValidationType.STRICT;

    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(doc="File to which the raw data will be written as JSON", fullName = "rawData", shortName = "rd", optional = true)
    public String jsonFile = null;

    @Argument(fullName = "additionalReportFile", shortName = "arf", doc="A text file listing the set of reports to display.  See Advanced Usage in the tool documentation", optional=true)
    public File additionalReportFile = null;

    private SampleDB sampleDB = null;

    protected ReportConfig[] getStandardWrappers() {
        PivotingTransformer transformer1 = new PivotingTransformer("CountVariants", Arrays.asList("Sample"), Arrays.asList(new PivotingTransformer.Pivot("FilterType", "nVariantLoci", null)));
        PivotingTransformer transformer2 = new PivotingTransformer("CountVariants", Arrays.asList("Sample"), Arrays.asList(new PivotingTransformer.Pivot("Contig", "nVariantLoci", null)), true);
        PivotingTransformer transformer3 = new PivotingTransformer("CountVariants", Arrays.asList("Contig"), Arrays.asList(new PivotingTransformer.Pivot("FilterType", "nVariantLoci", null)));
        PivotingTransformer transformer4 = new PivotingTransformer("CountVariants", Arrays.asList("EvalFeatureInput"), Arrays.asList(new PivotingTransformer.Pivot("FilterType", "nCalledLoci", null)));

        ReportConfig[] standardWrappers = new ReportConfig[]{
            new ReportConfig(new String[]{"EvalFeatureInput"}, TableReportDescriptor.getCountVariantsTable("Entire VCF", true)),
            new ReportConfig(new String[]{"EvalFeatureInput"}, BarPlotReportDescriptor.getVariantTypeBarPlot("Entire VCF")),
            new ReportConfig(new String[]{"EvalFeatureInput"}, TableReportDescriptor.getIndelTable("Entire VCF")),
            new ReportConfig(new String[]{"EvalFeatureInput"}, new TableReportDescriptor("Ti/Tv Data", "Entire VCF", "TiTvVariantEvaluator")),
            new ReportConfig(new String[]{"EvalFeatureInput"}, new TableReportDescriptor("Genotype Summary", "Entire VCF", "GenotypeFilterSummary")),

            new ReportConfig(new String[]{"Contig"}, TableReportDescriptor.getCountVariantsTable("By Contig", true)),
            new ReportConfig(new String[]{"Contig"}, BarPlotReportDescriptor.getVariantTypeBarPlot("By Contig")),
            new ReportConfig(new String[]{"Contig"}, TableReportDescriptor.getIndelTable("By Contig")),
            new ReportConfig(new String[]{"Contig"}, new TableReportDescriptor("Genotype Summary", "By Contig", "GenotypeFilterSummary", Arrays.asList("all"))),

            new ReportConfig(new String[]{"Sample"}, TableReportDescriptor.getCountVariantsTable("By Sample", true)),
            new ReportConfig(new String[]{"Sample"}, BarPlotReportDescriptor.getVariantTypeBarPlot("By Sample")),
            new ReportConfig(new String[]{"Sample"}, TableReportDescriptor.getIndelTable("By Sample")),
            new ReportConfig(new String[]{"Sample"}, new TableReportDescriptor("Ti/Tv Data", "By Sample", "TiTvVariantEvaluator", Arrays.asList("all"))),
            new ReportConfig(new String[]{"Sample"}, new TableReportDescriptor("Genotype Summary", "By Sample", "GenotypeFilterSummary", Arrays.asList("all"))),

            new ReportConfig(new String[]{"Sample", "FilterType"}, new TableReportDescriptor("Sites By Filter", "By Sample", "CountVariants", Arrays.asList("all"), transformer1)),
            new ReportConfig(new String[]{"Sample", "FilterType"}, BarPlotReportDescriptor.getSiteFilterTypeBarPlot("By Sample", transformer1)),

            new ReportConfig(new String[]{"Sample", "Contig"}, new TableReportDescriptor("Variants Per Contig", "By Sample", "CountVariants", Arrays.asList("all"), transformer2)),

            new ReportConfig(new String[]{"Contig", "FilterType"}, new TableReportDescriptor("Sites By Filter", "By Contig", "CountVariants", Arrays.asList("all"), transformer3)),
            new ReportConfig(new String[]{"Contig", "FilterType"}, BarPlotReportDescriptor.getSiteFilterTypeBarPlot("By Contig", transformer3)),

            new ReportConfig(new String[]{"FilterType"}, TableReportDescriptor.getCountVariantsTable("By Filter Type", true)),
            new ReportConfig(new String[]{"FilterType"}, BarPlotReportDescriptor.getVariantTypeBarPlot("By Filter Type")),

            new ReportConfig(new String[]{"EvalFeatureInput", "FilterType"}, new TableReportDescriptor("Variant Summary By Filter", "Entire VCF", "CountVariants", Arrays.asList("all"), transformer4)),
            new ReportConfig(new String[]{"EvalFeatureInput", "FilterType"}, BarPlotReportDescriptor.getSiteFilterTypeBarPlot("Entire VCF", transformer4))
        };

        return standardWrappers;
    }

    private Collection<VariantEvalWrapper> wrappers = new ArrayList<>();

    private Collection<VariantEvalWrapper> initializeReports()  {
        List<ReportConfig> configs = new ArrayList<>();
        configs.addAll(Arrays.asList(getStandardWrappers()));

        if (additionalReportFile != null) {
            IOUtil.assertFileIsReadable(additionalReportFile);

            configs.addAll(parseReportFile(additionalReportFile));
        }

        //preserve order to make resulting JSON consistent for test purposes
        Map<String, VariantEvalWrapper> reports = new LinkedHashMap<>();
        for (ReportConfig rc : configs) {
            String key = rc.getWrapperKey();
            VariantEvalWrapper wrapper = reports.getOrDefault(key, new VariantEvalWrapper(rc.stratifiers));
            wrapper.addReport(rc.rd);

            reports.put(key, wrapper);
        }

        return reports.values();
    }

    private List<ReportConfig> parseReportFile(File input) {
        List<ReportConfig> ret = new ArrayList<>();

        VCFHeader header = getHeaderForVariants();
        Map<String, Class<? extends VariantStratifier>> classMap = getStratifierClassMap();

        try (CSVReader reader = new CSVReader(IOUtil.openFileForBufferedUtf8Reading(input), '\t')) {
            String[] line;
            int i = 0;
            while ((line = reader.readNext()) != null) {
                i++;
                if (line[0].startsWith("#")) {
                    continue;
                }

                if (line.length < 4) {
                    throw new UserException.BadInput("Expected 4 columns on line: " + i + " of report config file");
                }

                String sectionLabel = StringUtils.trimToNull(line[0]);
                String reportLabel = StringUtils.trimToNull(line[1]);
                String stratifiers = StringUtils.trimToNull(line[2]);
                String infoField = StringUtils.trimToNull(line[3]);

                if (sectionLabel == null || reportLabel == null || stratifiers == null || infoField == null) {
                    throw new UserException.BadInput("Field " + infoField + " not found: on line " + i + " of report config file");
                }
                
                VCFInfoHeaderLine headerLine = header.getInfoHeaderLine(infoField);
                if (headerLine == null) {
                    throw new UserException.BadInput("INFO field " + infoField + " not found in VCF header (line " + i + " of report config file)");
                }

                if (headerLine.getType() != VCFHeaderLineType.Character && headerLine.getType() != VCFHeaderLineType.String && headerLine.getType() != VCFHeaderLineType.Integer) {
                    throw new UserException.BadInput("Field " + infoField + " was not a supported type (line " + i + " of report config file).  Currently supported types are Character, String and Integer");
                }

                TableReportDescriptor rd = new TableReportDescriptor.InfoFieldTableReportDescriptor(reportLabel, sectionLabel, infoField);
                String[] stratList = stratifiers.split(",");
                for (String strat : stratList) {
                    if (!classMap.containsKey(strat)) {
                        throw new UserException.BadInput("Stratifier not found " + strat + " (line " + i + " of report config file)");
                    }
                }
                
                ret.add(new ReportConfig(stratList, rd));
            }
        }
        catch (IOException e) {
            throw new UserException.BadInput("Unable to parse report file", e);
        }

        return ret;
    }

    public class ReportConfig {
        ArrayList<String> stratifiers;
        ReportDescriptor rd;

        public ReportConfig(String[] stratifiers, ReportDescriptor rd) {
            this.stratifiers = new ArrayList<>(Arrays.asList(stratifiers));
            this.rd = rd;
        }

        public String getWrapperKey() {
            return StringUtils.join(this.stratifiers, ";");
        }
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        Utils.nonNull(outFile);

        if (jsonFile != null) {
            File json = new File(jsonFile);
            IOUtil.assertFileIsWritable(json);
        }

        sampleDB = initializeSampleDB();

        this.wrappers = initializeReports();

        //configure the child walkers
        for (VariantEvalWrapper wrapper : this.wrappers){
            wrapper.configureWalker(this);

            wrapper.walker.onTraversalStart();
        }
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        for (VariantEvalWrapper wrapper : this.wrappers) {
            wrapper.walker.apply(variant, readsContext, referenceContext, featureContext);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        for (VariantEvalWrapper wrapper : this.wrappers){
            wrapper.walker.onTraversalSuccess();
        }

        Map<String, SectionJsonDescriptor> sectionMap = new LinkedHashMap<>();
        Map<String, Class<? extends VariantStratifier>> classMap = getStratifierClassMap();

        for (VariantEvalWrapper wrapper : this.wrappers) {
            try (BufferedReader sampleReader = IOUtil.openFileForBufferedUtf8Reading(wrapper.outFile)) {
                sampleReader.readLine(); //read first GATKReport line

                for (int i=0; i<wrapper.evaluationModules.size(); i++){
                    //NOTE: this output will have one table per eval module. Iterate
                    GATKReportTable table = new GATKReportTable(sampleReader, GATKReportVersion.V1_1);
                    List<ReportDescriptor> rds = wrapper.getReportsForModule(table.getTableName());
                    Map<String, String> descriptionMap = new HashMap<>();
                    Class<? extends VariantStratifier> evalClass = classMap.get(table.getTableName());
                    if (evalClass != null){
                        AnalysisModuleScanner scanner = new AnalysisModuleScanner(evalClass);
                        Map<Field, DataPoint> fieldDataPointMap = scanner.getData();
                        for (Field f : fieldDataPointMap.keySet()){
                            descriptionMap.put(f.getName(), fieldDataPointMap.get(f).description());
                        }
                    }
                    if (rds.isEmpty()){
                        throw new GATKException("No report registered for GATK table: " + table.getTableName());
                    }

                    for (ReportDescriptor rd : rds){
                        if (!sectionMap.containsKey(rd.sectionLabel)){
                            sectionMap.put(rd.sectionLabel, new SectionJsonDescriptor(rd.sectionLabel, wrapper.stratifications));
                        }

                        sectionMap.get(rd.sectionLabel).addReportDescriptor(rd, table, descriptionMap, sampleDB);
                    }
                }
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            wrapper.outFile.delete();
        }

        try {
            List<SectionJsonDescriptor> sections = new ArrayList<>();

            sectionMap.keySet().forEach(key -> sections.add(sectionMap.get(key)));

            try (PrintWriter writer = new PrintWriter(IOUtil.openFileForBufferedWriting(new File(outFile))); PrintWriter jsonWriter = (jsonFile == null ? null : new PrintWriter(IOUtil.openFileForBufferedWriting(new File(jsonFile))))) {
                HtmlGenerator generator = new HtmlGenerator();
                generator.generateHtml(sections, writer, jsonWriter);
            }
        }
        catch (IOException e){
            throw new GATKException(e.getMessage(), e);
        }

        return super.onTraversalSuccess();
    }

    public static class VariantEvalWrapper {
        private VariantEvalChild walker;
        //private ByteArrayOutputStream out = new ByteArrayOutputStream();
        private File outFile = IOUtils.createTempFile("variantQC_data", ".txt");
        List<String> stratifications;
        Set<String> evaluationModules = new HashSet<>();
        List<String> infoFieldEvaluators = new ArrayList<>();
        List<ReportDescriptor> reportDescriptors = new ArrayList<>();

        public VariantEvalWrapper(List<String> stratifications) {
            this.stratifications = new ArrayList<>();
            this.stratifications.add("Filter");
            this.stratifications.addAll(stratifications);
        }

        public void addReport(ReportDescriptor rd) {
            this.reportDescriptors.add(rd);

            if (rd instanceof TableReportDescriptor.InfoFieldTableReportDescriptor) {
                infoFieldEvaluators.add(((TableReportDescriptor.InfoFieldTableReportDescriptor)rd).getInfoFieldName());
            }
            else {
                this.evaluationModules.add(rd.getEvaluatorModuleName());
            }
        }

        public List<ReportDescriptor> getReportsForModule(String evalModule){
            List<ReportDescriptor> ret = new ArrayList<>();
            for (ReportDescriptor rd : reportDescriptors){
                if (rd.evaluatorModuleName.equals(evalModule)){
                    ret.add(rd);
                }
            }

            return ret;
        }

        public File getOutFile() {
            return outFile;
        }

        public void configureWalker(VariantQC variantQC) {
            this.walker = new VariantEvalChild(variantQC, this, variantQC.getDrivingVariantsFeatureInput(), infoFieldEvaluators);
        }
    }

    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(pedigreeValidationType);
        if (pedigreeFile != null)
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));

        Collection<String> samples = getHeaderForVariants().getSampleNamesInOrder();
        if (samples != null) {
            sampleDBBuilder.addSamplesFromSampleNames(samples);
        }

        return sampleDBBuilder.getFinalSampleDB();
    }

    //TODO: add getter to VariantEvalUtils
    private Map<String, Class<? extends VariantStratifier>> getStratifierClassMap() {
        Map<String, Class<? extends VariantStratifier>> stratifierClasses = new HashMap<>();
        Reflections reflectionsStrat = new Reflections("org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications");
        Set<Class<? extends VariantStratifier>> allClasses = reflectionsStrat.getSubTypesOf(VariantStratifier.class);
        for (Class<? extends VariantStratifier> clazz : allClasses) {
            stratifierClasses.put(clazz.getSimpleName(), clazz);
        }

        return stratifierClasses;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }
}
