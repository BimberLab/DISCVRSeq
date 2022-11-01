package com.github.discvrseq.walkers.variantqc;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvValidationException;
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
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.AnalysisModuleScanner;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.report.GATKReportVersion;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.stream.Collectors;

/**
 * VariantQC will generate a user-friendly HTML report that aggregates data from a VCF file by site, filter status, sample, and other stratifiers in order
 * to provide various summary tables and graphs.  In most cases the tables overlay bar graphs for numeric columns in order to help identify outliers.
 * This tool is analogous to FASTQC or MultiQC, which provide similar HTML summary reports for different types of sequence data.
 *
 * <h3>Usage examples:</h3>
 *
 *  Note: The public GATK Reference Bundle can be downloaded to provide example VCFs and a pre-indexed GRCh38 human genome build.  <a href="https://software.broadinstitute.org/gatk/documentation/article.php?id=11017">See here for download instructions.</a>
 *
 *  We have also provide a small simplified VCF.  You can <a href="SimpleExample.vcf.gz">download the VCF</a> and <a href="SimpleExample.vcf.gz.tbi">download the VCF index</a>, with can be used with the <a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/">human_g1k_v37 genome</a>.
 *
 * <h4>Basic Usage, Without Pedigree:</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R human_g1k_v37.fasta \
 *     -V SimpleExample.vcf.gz \
 *     -O output.html
 * </pre>
 *
 * <h4>Report With Pedigree (required for gender information to display):</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R human_g1k_v37.fasta \
 *     -ped myPedigree.ped \
 *     -V input.vcf \
 *     -O output.html
 * </pre>
 *
 * <h4>Report With Multiple Input VCFs (not well tested). Note: you must supply a name for each VCF in the argument:</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R human_g1k_v37.fasta \
 *     -ped myPedigree.ped \
 *     -V:vcf1 input1.vcf \
 *     -V:vcf2 input2.vcf \
 *     -O output.html
 * </pre>
 *
 * <h3>Variant QC Output Report:</h3>
 * Below is an image of an example report, showing a VCF summarized by sample.  The left-hand navigation allows the user to toggle between different stratifications (entire VCF, sample, contig, etc.).
 * A complete example HTML report can be <a href="resources/variantQCSampleReport.html">viewed here</a>.  The report contains a series of tables or graphs summarizing different aspects of the data.  Below are the main types:
 * <br>
*  <ul>
 *      <li>Variant Summary: A table displaying a summary of the total variants by type (SNP, insertion, deletion, MNP, etc.), </li>
 *      <li>Variant Type: a bar plot summarizing variant by type (SNP, insertion, deletion, MNP, etc.)</li>
 *      <li>Genotype Summary: A table summarizing the total called/non-called genotypes</li>
 *      <li>SNP/Indel Summary:</li>
 *      <li>Ti/Tv Data: A table displaying a summary of <a href="https://en.wikipedia.org/wiki/Transversion">transition and transversion</a> mutations in the dataset</li>
 *      <li>Filter Type: A bar plot summarizing variants by filter</li>
*  </ul>
 * <p></p>
 * <a href="resources/variantQCSampleReport.html"><img src="images/variantQC_ExampleReport.jpg" style="width: 75%" title="Click to view a sample report"/></a>
 *
 * <h3>Advanced Usage:</h3>
 *
 * VariantQC is designed to produce general-purpose reports applicable to a range of users; however, several mechanisms may make the output more useful to your data:
 *
 * <h4>Write the raw output data to a separate file as JSON, which can be parsed separately:</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R human_g1k_v37.fasta \
 *     -V input.vcf \
 *     -rd output.json \
 *     -O output.html
 * </pre>
 *
 * <h4>Define additional reports to be included in the HTML report:</h4>
 * <pre>
 * java -jar DISCVRSeq.jar VariantQC \
 *     -R human_g1k_v37.fasta \
 *     -V input.vcf \
 *     --additionalReportFile reports.txt \
 *     -O output.html
 * </pre>
 *
 * Where reports.txt is a tab-delimited file with one report per line and no header.  Comment lines can be included, beginning with '#'.  This file has 4 columns:
 * <br>
 * <ul>
 * <li>Section label: Corresponds to the report group.  This is typically 'Entire VCF', 'By Contig', etc.; however, any value can be used.</li>
 * <li>Report label: The title used for this report.  For example: 'Variant Summary'</li>
 * <li>Stratifications: A comma-separated list of the stratifications to use when aggregating data.  Allowable values are: VCF, Sample, Contig, and Filter</li>
 * <li>INFO field name: The name of the INFO field attribute to aggregate.  This field must be of type character, string or integer.  Note: this will produce a table summarizing the total variants for each level of this variable.  Therefore a field with a wide range of possible values may not be appropriate to summarize in this manner.  For example, while integer fields are support since in certain cases the value is bounded and will have a reasonable number of unique values.</li>
 * </ul>
 * <br>
 * Example Report File:
 * <pre>
 * By Sample	Example Report	Sample	PURPOSE
 * By Contig	Example Report2	Sample,Contig	PURPOSE
 * </pre>
 *
 * <h4>Other Usage Suggestions:</h4>
 * Upstream processing of your VCF can enhance the value of the VariantQC report for your data. Our group routinely performs quality filtering on our VCFs, which saves information about the filter type in the FILTER field (<a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php"></a>see VariantFiltration</a>).  FilterType is used to stratify data in VariantQC, allowing us to view sample, VCF, or chromosome differences.
 * <br><br>
 * While VariantQC does not direct support automatic filtering/flagging of VCFs or samples based on thresholds, our group accomplishes this by writing a separate script to read the raw JSON output (see --rawData argument), and then applying our application-specific criteria.
 *
 * <h3>Authors:</h3>
 * <ul>
 *  <li>Ben Bimber, Ph.D.</li>
 *  <li>Melissa Yan</li>
 * </ul>
 *
 * <h3>Acknowledgements:</h3>
 * This tool internally uses <a href="https://github.com/broadinstitute/gatk">GATK4's VariantEval</a> to aggregate data, and borrows heavily from <a href="https://multiqc.info/">MultiQC</a> for the HTML in the final report.
 * We thank Philip Ewels for the use of MultiQC HTML/CSS/JS templates, and Chris Norman of the Broad Institute for assistance with VariantEval.
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This will generate an HTML summary report for a given VCF, analogous to the report FastQC creates for raw read data.",
        oneLineSummary = "Generate HTML summary report for a VCF",
        programGroup = VariantManipulationProgramGroup.class
)
public class VariantQC extends MultiVariantWalkerGroupedOnStart {

    @Argument(fullName = StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName = StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for identifying Mendelian violations", optional=true)
    private GATKPath pedigreeFile;

    @Argument(fullName = "pedigreeValidationType", shortName = "pedValidationType", doc="The strictness for validating the pedigree.  Can be either STRICT or SILENT.  Default is STRICT", optional=true)
    private PedigreeValidationType pedigreeValidationType = PedigreeValidationType.STRICT;

    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(doc="File to which the raw data will be written as JSON", fullName = "rawData", shortName = "rd", optional = true)
    public String jsonFile = null;

    @Argument(fullName = "additionalReportFile", shortName = "arf", doc="This is an advanced usage feature.  The user can define additional reports to display.  Each report will read the value of the supplied INFO field, and make a table summarizing the number of variants for each value of this field.  For example, if your VCF has the INFO field '', and this has values of , a table will be created summarizing the number of variants for each state.  These will be stratified as defined in your file.  Only INFO fields of type character, string and integer can be used.  The file itself should be a tab-delimited file with one report per line, no header, and 4 columns: Section Label, Report Label, Stratification(s), and name of the INFO field to summarize.  The TSV file is described in greater detail in the Advanced Usage section.", optional=true)
    public File additionalReportFile = null;

    @Argument(fullName = "maxContigs", shortName = "maxContigs", doc="Many VariantQC reports stratify data by contig.  If the genome contains a large number of chromosomes, such as lots of unplaced contigs, this can bog down these reports in the final HTML file. As a default, VariantQC will only process the first 40 contigs, by length. This can be increased using this argument. See also --contigsToRetain", optional=true)
    public int maxContigs = 40;

    @Argument(fullName = "contigsToRetain", doc="If --maxContigs is used, the first X contigs, are retained, sorted by length and preferentially retaining the longest contigs. This can be used to specify one or more additional contigs that are retained, even if they would otherwise be removed.", optional=true)
    public List<String> contigsToRetain = new ArrayList<>(Collections.singletonList("MT"));

    @Argument(fullName = "threads", doc="The number of threads to use.", optional=true)
    public int threads = 1;

    private SampleDB sampleDB = null;

    private List<ReportConfig> getStandardWrappers(boolean hasSamples, boolean isMultiVcf) {
        PivotingTransformer transformer1 = new PivotingTransformer("CountVariants", Arrays.asList("Sample"), isMultiVcf, Arrays.asList(new PivotingTransformer.Pivot("FilterType", "nVariantLoci", null)));
        PivotingTransformer transformer2 = new PivotingTransformer("CountVariants", Arrays.asList("Sample"), isMultiVcf, Arrays.asList(new PivotingTransformer.Pivot("Contig", "nVariantLoci", null)), true);
        PivotingTransformer transformer3 = new PivotingTransformer("CountVariants", Arrays.asList("Contig"), isMultiVcf, Arrays.asList(new PivotingTransformer.Pivot("FilterType", "nVariantLoci", null)));
        PivotingTransformer transformer4 = new PivotingTransformer("CountVariants", Arrays.asList("EvalFeatureInput"), isMultiVcf, Arrays.asList(new PivotingTransformer.Pivot("FilterType", "nCalledLoci", null)));

        List<ReportConfig> standardWrappers = new ArrayList<>(Arrays.asList(
                new ReportConfig(Arrays.asList("EvalFeatureInput"), TableReportDescriptor.getCountVariantsTable("Entire VCF", isMultiVcf, true)),
                new ReportConfig(Arrays.asList("EvalFeatureInput"), BarPlotReportDescriptor.getVariantTypeBarPlot("Entire VCF", isMultiVcf)),
                new ReportConfig(Arrays.asList("EvalFeatureInput"), TableReportDescriptor.getIndelTable("Entire VCF", isMultiVcf)),
                new ReportConfig(Arrays.asList("EvalFeatureInput"), new TableReportDescriptor("Ti/Tv Data", "Entire VCF", isMultiVcf, "TiTvVariantEvaluator")),
                new ReportConfig(Arrays.asList("EvalFeatureInput"), new TableReportDescriptor("Genotype Summary", "Entire VCF", isMultiVcf, "GenotypeFilterSummary")),

                new ReportConfig(Arrays.asList("Contig"), TableReportDescriptor.getCountVariantsTable("By Contig", isMultiVcf, true)),
                new ReportConfig(Arrays.asList("Contig"), BarPlotReportDescriptor.getVariantTypeBarPlot("By Contig", isMultiVcf)),
                new ReportConfig(Arrays.asList("Contig"), TableReportDescriptor.getIndelTable("By Contig", isMultiVcf)),
                new ReportConfig(Arrays.asList("Contig"), new TableReportDescriptor("Genotype Summary", "By Contig", isMultiVcf, "GenotypeFilterSummary", Arrays.asList("all")))
        ));

        if (hasSamples) {
            standardWrappers.add(new ReportConfig(Arrays.asList("Sample"), TableReportDescriptor.getCountVariantsTable("By Sample", isMultiVcf, true)));
            standardWrappers.add(new ReportConfig(Arrays.asList("Sample"), BarPlotReportDescriptor.getVariantTypeBarPlot("By Sample", isMultiVcf)));
            standardWrappers.add(new ReportConfig(Arrays.asList("Sample"), TableReportDescriptor.getIndelTable("By Sample", isMultiVcf)));
            standardWrappers.add(new ReportConfig(Arrays.asList("Sample"), new TableReportDescriptor("Ti/Tv Data", "By Sample", isMultiVcf, "TiTvVariantEvaluator", Arrays.asList("all"))));
            standardWrappers.add(new ReportConfig(Arrays.asList("Sample"), new TableReportDescriptor("Genotype Summary", "By Sample", isMultiVcf, "GenotypeFilterSummary", Arrays.asList("all"))));

            standardWrappers.add(new ReportConfig(Arrays.asList("Sample", "FilterType"), new TableReportDescriptor("Sites By Filter", "By Sample", isMultiVcf, "CountVariants", Arrays.asList("all"), transformer1)));
            standardWrappers.add(new ReportConfig(Arrays.asList("Sample", "FilterType"), BarPlotReportDescriptor.getSiteFilterTypeBarPlot("By Sample", isMultiVcf, transformer1)));

            standardWrappers.add(new ReportConfig(Arrays.asList("Sample", "Contig"), new TableReportDescriptor("Variants Per Contig", "By Sample", isMultiVcf, "CountVariants", Arrays.asList("all"), transformer2)));
        }

        standardWrappers.add(new ReportConfig(Arrays.asList("Contig", "FilterType"), new TableReportDescriptor("Sites By Filter", "By Contig", isMultiVcf, "CountVariants", Arrays.asList("all"), transformer3)));
        standardWrappers.add(new ReportConfig(Arrays.asList("Contig", "FilterType"), BarPlotReportDescriptor.getSiteFilterTypeBarPlot("By Contig", isMultiVcf, transformer3)));

        standardWrappers.add(new ReportConfig(Arrays.asList("FilterType"), TableReportDescriptor.getCountVariantsTable("By Filter Type", isMultiVcf, true)));
        standardWrappers.add(new ReportConfig(Arrays.asList("FilterType"), BarPlotReportDescriptor.getVariantTypeBarPlot("By Filter Type", isMultiVcf)));

        standardWrappers.add(new ReportConfig(Arrays.asList("EvalFeatureInput", "FilterType"), new TableReportDescriptor("Variant Summary By Filter", "Entire VCF", isMultiVcf, "CountVariants", Arrays.asList("all"), transformer4)));
        standardWrappers.add(new ReportConfig(Arrays.asList("EvalFeatureInput", "FilterType"), BarPlotReportDescriptor.getSiteFilterTypeBarPlot("Entire VCF", isMultiVcf, transformer4)));

        return standardWrappers;
    }

    private Collection<VariantEvalWrapper> wrappers = new ArrayList<>();

    private Collection<VariantEvalWrapper> initializeReports()  {
        List<ReportConfig> configs = new ArrayList<>();
        configs.addAll(getStandardWrappers(!getSamplesForVariants().isEmpty(), isMultiVcf()));

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

    private boolean isMultiVcf() {
        return getDrivingVariantsFeatureInputs().size() > 1;
    }

    private List<ReportConfig> parseReportFile(File input) {
        List<ReportConfig> ret = new ArrayList<>();

        VCFHeader header = getHeaderForVariants();
        Map<String, Class<? extends VariantStratifier>> classMap = new HashMap<>(VariantEvalEngine.getStratifierClasses());
        classMap.put(Contig.class.getSimpleName(), Contig.class);

        try (CSVReader reader = CsvUtils.getTsvReader(input)) {
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

                TableReportDescriptor rd = new TableReportDescriptor.InfoFieldTableReportDescriptor(reportLabel, sectionLabel, isMultiVcf(), infoField);
                List<String> stratList = new ArrayList<>(Arrays.asList(stratifiers.split(",")));

                //allow user-friendly translation:
                ListIterator<String> it = stratList.listIterator();
                while (it.hasNext()) {
                    String v = it.next();
                    if ("VCF".equals(v)) {
                        it.set("EvalFeatureInput");
                    }
                }

                for (String strat : stratList) {

                    if (!classMap.containsKey(strat)) {
                        Set<String> allowable = new TreeSet<>(classMap.keySet());
                        allowable.add("VCF");
                        allowable.remove("EvalFeatureInput");

                        throw new UserException.BadInput("Stratifier not found " + strat + " (line " + i + " of report config file).  Allowable values are: " + StringUtils.join(allowable, ", "));
                    }
                }

                ret.add(new ReportConfig(stratList, rd));
            }
        }
        catch (CsvValidationException | IOException e) {
            throw new UserException.BadInput("Unable to parse report file", e);
        }

        return ret;
    }

    public static class ReportConfig {
        ArrayList<String> stratifiers;
        ReportDescriptor rd;

        public ReportConfig(List<String> stratifiers, ReportDescriptor rd) {
            this.stratifiers = new ArrayList<>(stratifiers);
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

        if (isMultiVcf()) {
            Set<String> unique = new HashSet<>();
            Set<String> duplicates = new HashSet<>();
            getDrivingVariantsFeatureInputs().stream().map(x -> x.hasUserSuppliedName() ? x.getName() : "eval").forEach(x -> {
                if (!unique.add(x)) {
                    duplicates.add(x);
                }
            });

            if (!duplicates.isEmpty()) {
                throw new UserException("When supplying multiple input VCFs, they must each be assigned a unique name on the command line, such as: '--variant:vcf1 myVcf.vcf.gz'");
            }
        }

        this.wrappers = initializeReports();
        logger.info("Total VariantEval instances: " + this.wrappers.size());

        //configure the child walkers
        this.getTraversalIntervals(); //initialize potential contig override
        for (VariantEvalWrapper wrapper : this.wrappers){
            wrapper.configureEngine(this);
        }

        if (threads > 1) {
            executor = Executors.newScheduledThreadPool(threads);
        }
    }

    ScheduledExecutorService executor = null;

    @Override
    public void apply(final List<VariantContext> list, final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
        if (executor != null) {
            List<ApplyRunner> toRun = new ArrayList<>();
            for (final VariantEvalWrapper wrapper : this.wrappers) {
                toRun.add(new ApplyRunner(list, referenceContext, wrapper));
            }

            try {
                executor.invokeAll(toRun);
            }
            catch (InterruptedException e) {
                throw new IllegalStateException("Error running VariantQC", e);
            }
        }
        else {
            for (VariantEvalWrapper wrapper : this.wrappers) {
                wrapper.engine.apply(list, referenceContext);
            }
        }
    }

    public static class ApplyRunner implements Callable<Boolean> {
        final List<VariantContext> list;
        final ReferenceContext referenceContext;
        final VariantEvalWrapper wrapper;

        public ApplyRunner(final List<VariantContext> list, final ReferenceContext referenceContext, final VariantEvalWrapper wrapper) {
            // NOTE: ensureAnnotations() mutates the VC, so call this one time, upfront:
            this.list = list.stream().map(vc -> wrapper.engine.ensureAnnotations(vc, vc)).collect(Collectors.toList());
            this.referenceContext = new ReferenceContext(referenceContext, referenceContext.getInterval());
            this.wrapper = wrapper;
        }

        @Override
        public Boolean call() throws Exception {
            wrapper.engine.apply(list, referenceContext);

            return true;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (executor != null) {
            executor.shutdown();
        }

        //TODO: option to write to disk
        for (VariantEvalWrapper wrapper : this.wrappers){
            wrapper.engine.finalizeReport(wrapper.getOutFile());
        }

        Map<String, SectionJsonDescriptor> sectionMap = new LinkedHashMap<>();
        Map<String, Class<? extends VariantEvaluator>> classMap = VariantEvalEngine.getEvaluatorClasses();

        for (VariantEvalWrapper wrapper : this.wrappers) {
            try (BufferedReader sampleReader = IOUtil.openFileForBufferedUtf8Reading(wrapper.outFile)) {
                sampleReader.readLine(); //read first GATKReport line

                for (int i=0; i<wrapper.getEvaluationModules().size(); i++){
                    //NOTE: this output will have one table per eval module. Iterate
                    GATKReportTable table = new GATKReportTable(sampleReader, GATKReportVersion.V1_1);
                    List<ReportDescriptor> rds = wrapper.getReportsForModule(table.getTableName());
                    Map<String, String> descriptionMap = new HashMap<>();
                    Class<? extends VariantEvaluator> evalClass = classMap.get(table.getTableName());

                    //TODO improve this when refactoring VariantEvalEngine
                    if (evalClass == null && table.getTableName().startsWith(InfoFieldEvaluator.class.getSimpleName() + "-")) {
                        evalClass = InfoFieldEvaluator.class;
                    }

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

            try (PrintWriter writer = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(new File(outFile))); PrintWriter jsonWriter = (jsonFile == null ? null : new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(new File(jsonFile))))) {
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
        private VariantEvalEngine engine;
        //private ByteArrayOutputStream out = new ByteArrayOutputStream();
        private File outFile = IOUtils.createTempFile("variantQC_data", ".txt");

        private List<String> stratifications;
        private Set<String> evaluationModules = new HashSet<>();
        private List<String> infoFields = new ArrayList<>();
        private List<ReportDescriptor> reportDescriptors = new ArrayList<>();

        public VariantEvalWrapper(List<String> stratifications) {
            this.stratifications = new ArrayList<>();
            this.stratifications.add("Filter");
            this.stratifications.addAll(stratifications);
        }

        public void addReport(ReportDescriptor rd) {
            this.reportDescriptors.add(rd);

            if (rd instanceof TableReportDescriptor.InfoFieldTableReportDescriptor) {
                infoFields.add(((TableReportDescriptor.InfoFieldTableReportDescriptor)rd).getInfoFieldName());
            }
            else {
                this.evaluationModules.add(rd.getEvaluatorModuleName());
            }
        }

        public Set<String> getEvaluationModules() {
            Set<String> ret = new TreeSet<>(evaluationModules);
            infoFields.forEach(x -> ret.add(TableReportDescriptor.InfoFieldTableReportDescriptor.getEvalModuleSimpleName(x)));

            return ret;
        }

        public List<ReportDescriptor> getReportsForModule(String evalModule){
            List<ReportDescriptor> ret = new ArrayList<>();
            for (ReportDescriptor rd : reportDescriptors){
                if (rd.getEvaluatorModuleName().equals(evalModule)){
                    ret.add(rd);
                }
            }

            return ret;
        }

        public File getOutFile() {
            return outFile;
        }

        public void configureEngine(VariantQC variantQC) {

            VariantEvalArgumentCollection args = new VariantEvalArgumentCollection();
            args.evals = Collections.unmodifiableList(variantQC.getDrivingVariantsFeatureInputs());
            args.modulesToUse = new ArrayList<>(this.evaluationModules);
            args.noStandardModules = true;
            args.stratificationsToUse = stratifications;
            args.noStandardStratifications = true;
            args.pedigreeFile = variantQC.pedigreeFile;
            args.pedigreeValidationType = variantQC.pedigreeValidationType;

            this.engine = new ExtendedVariantEvalEngine(args, variantQC.features, variantQC.getTraversalIntervals(), variantQC.getSequenceDictionaryForDrivingVariants(), variantQC.getSamplesForVariants(), infoFields, variantQC.maxContigs, variantQC.contigsToRetain);
        }
    }

    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(pedigreeValidationType);
        if (pedigreeFile != null)
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));

        Collection<String> samples = getSamplesForVariants();
        if (samples != null) {
            sampleDBBuilder.addSamplesFromSampleNames(samples);
        }

        return sampleDBBuilder.getFinalSampleDB();
    }

    @Override
    public boolean requiresReference() {
        return true;
    }
}
