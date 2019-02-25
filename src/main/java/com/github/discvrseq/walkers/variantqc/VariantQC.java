package com.github.discvrseq.walkers.variantqc;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.AnalysisModuleScanner;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.utils.SimpleInterval;
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
 * contains a series of tables or graphs summarizing different aspects of the data.
 * <p></p>
 * <img src="images/variantQCSampleReport1.png" style="width: 75%"/>
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

    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    private PrintWriter writer;

    private SampleDB sampleDB = null;

    protected VariantEvalWrapper[] wrappers = new VariantEvalWrapper[]{
            // TODO: we need to add evaluators that track specific JEXL expressions, such as:
            // # mendelian violations / site, which we plot as a histogram.  perhaps also AF, as histogram
            // We should be able to implement this as either:
            // a) some type of generic VariantEvaluator class accepting a JEXL expression that will return a number.  the problem here is we dont know the states upfront.
            // b) an alternative could be to skip VariantEval and write a new standalone HistogramWalker class.  We would need to maintain a separate internal set of these; however, I could imagine this tool being generically useful.

            new VariantEvalWrapper("Entire VCF", new String[]{"EvalFeatureInput"}, new String[]{"CountVariants", "IndelSummary", "TiTvVariantEvaluator", "GenotypeFilterSummary"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(true),
                    BarPlotReportDescriptor.getVariantTypeBarPlot(),
                    TableReportDescriptor.getIndelTable(),
                    new TableReportDescriptor("Ti/Tv Data", "TiTvVariantEvaluator"),
                    new TableReportDescriptor("Genotype Summary", "GenotypeFilterSummary")
            }),
            new VariantEvalWrapper("By Contig", new String[]{"Contig"}, new String[]{"CountVariants", "IndelSummary", "GenotypeFilterSummary"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(true),
                    BarPlotReportDescriptor.getVariantTypeBarPlot(),
                    TableReportDescriptor.getIndelTable(),
                    new TableReportDescriptor("Genotype Summary", "GenotypeFilterSummary", Arrays.asList("all"))
            }),
            new VariantEvalWrapper("By Sample", new String[]{"Sample"}, new String[]{"CountVariants", "IndelSummary", "TiTvVariantEvaluator", "GenotypeFilterSummary"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(true),
                    BarPlotReportDescriptor.getVariantTypeBarPlot(),
                    TableReportDescriptor.getIndelTable(),
                    new TableReportDescriptor("Ti/Tv Data", "TiTvVariantEvaluator", Arrays.asList("all")),
                    new TableReportDescriptor("Genotype Summary", "GenotypeFilterSummary", Arrays.asList("all"))
            }),
            new VariantEvalWrapper("By Sample", new String[]{"Sample", "FilterType"}, new String[]{"CountVariants"}, new ReportDescriptor[]{
                    new TableReportDescriptor("Sites By Filter", "CountVariants", Arrays.asList("all")),
                    BarPlotReportDescriptor.getSiteFilterTypeBarPlot(),
            }, Arrays.asList(new PivotingTransformer("CountVariants", Arrays.asList("Sample"), Arrays.asList(
                    new PivotingTransformer.Pivot("FilterType", "nVariantLoci", null)
            )))),
            new VariantEvalWrapper("By Sample", new String[]{"Sample", "Contig"}, new String[]{"CountVariants"}, new ReportDescriptor[]{
                    new TableReportDescriptor("Variants Per Contig", "CountVariants", Arrays.asList("all")),
            }, Arrays.asList(new PivotingTransformer("CountVariants", Arrays.asList("Sample"), Arrays.asList(
                    new PivotingTransformer.Pivot("Contig", "nVariantLoci", null)
            ), true))),
            new VariantEvalWrapper("By Contig", new String[]{"Contig", "FilterType"}, new String[]{"CountVariants"}, new ReportDescriptor[]{
                    new TableReportDescriptor("Sites By Filter", "CountVariants", Arrays.asList("all")),
                    BarPlotReportDescriptor.getSiteFilterTypeBarPlot(),
            }, Arrays.asList(new PivotingTransformer("CountVariants", Arrays.asList("Contig"), Arrays.asList(
                    new PivotingTransformer.Pivot("FilterType", "nVariantLoci", null)
            )))),
            new VariantEvalWrapper("By Filter Type", new String[]{"FilterType"}, new String[]{"CountVariants"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(true),
                    BarPlotReportDescriptor.getVariantTypeBarPlot()
            }),
            new VariantEvalWrapper("Entire VCF", new String[]{"EvalFeatureInput", "FilterType"}, new String[]{"CountVariants"}, new ReportDescriptor[]{
                    new TableReportDescriptor("Variant Summary By Filter", "CountVariants", Arrays.asList("all")),
                    BarPlotReportDescriptor.getSiteFilterTypeBarPlot()
            }, Arrays.asList(new PivotingTransformer("CountVariants", Arrays.asList("EvalFeatureInput"), Arrays.asList(
                    new PivotingTransformer.Pivot("FilterType", "nCalledLoci", null)
            ))))
    };

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        Utils.nonNull(outFile);
        writer = new PrintWriter(IOUtil.openFileForBufferedWriting(new File(outFile)));

        sampleDB = initializeSampleDB();

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
        Map<String, Class<? extends VariantStratifier>> classMap = getEvaluationClassMap();

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

                    GATKReportTableTransformer transformer = wrapper.transformerMap.get(table.getTableName());
                    if (transformer != null){
                        table = transformer.transform(table, sampleDB);
                    }

                    if (!sectionMap.containsKey(wrapper.sectionLabel)){
                        sectionMap.put(wrapper.sectionLabel, new SectionJsonDescriptor(wrapper.sectionLabel, wrapper.stratifications));
                    }

                    for (ReportDescriptor rd : rds){
                        sectionMap.get(wrapper.sectionLabel).addReportDescriptor(rd, table, descriptionMap);
                    }
                }
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }
        }

        try {
            List<SectionJsonDescriptor> sections = new ArrayList<>();

            sectionMap.keySet().forEach(key -> sections.add(sectionMap.get(key)));

            HtmlGenerator generator = new HtmlGenerator();
            generator.generateHtml(sections, writer);
        }
        catch (IOException e){
            throw new GATKException(e.getMessage(), e);
        }

        writer.close();

        return super.onTraversalSuccess();
    }


    private static class VariantEvalChild extends VariantEval {
        private final VariantQC variantQC;

        public VariantEvalChild(VariantQC variantQC, VariantEvalWrapper wrapper){
            this.variantQC = variantQC;
            try
            {
                this.evals = Collections.singletonList(variantQC.getDrivingVariantsFeatureInput());
                this.outFile = wrapper.outFile;

                this.MODULES_TO_USE = wrapper.evaluationModules;
                this.NO_STANDARD_MODULES = true;
                this.STRATIFICATIONS_TO_USE = wrapper.stratifications;
                this.NO_STANDARD_STRATIFICATIONS = true;

                //TODO: set reference??

                this.onStartup();
            }
            catch (Exception e)
            {
                throw new GATKException(e.getMessage(), e);
            }
        }

        @Override
        public List<SimpleInterval> getTraversalIntervals() {
            return variantQC.getTraversalIntervals();
        }
    }

    private static class VariantEvalWrapper {
        private VariantEvalChild walker;
        //private ByteArrayOutputStream out = new ByteArrayOutputStream();
        private File outFile = IOUtils.createTempFile("variantQC_data", ".txt");
        List<String> stratifications;
        List<String> evaluationModules;
        String sectionLabel;
        ReportDescriptor[] reportDescriptors;
        Map<String, GATKReportTableTransformer> transformerMap;

        public VariantEvalWrapper(String sectionLabel, String[] stratifications, String[] evaluationModules, ReportDescriptor[] reportDescriptors){
            this(sectionLabel, stratifications, evaluationModules, reportDescriptors, null);
        }

        public VariantEvalWrapper(String sectionLabel, String[] stratifications, String[] evaluationModules, ReportDescriptor[] reportDescriptors, List<GATKReportTableTransformer> transformers) {
            this.stratifications = new ArrayList<>();
            this.stratifications.add("Filter");
            this.stratifications.addAll(Arrays.asList(stratifications));

            this.evaluationModules = Arrays.asList(evaluationModules);
            this.sectionLabel = sectionLabel;

            this.reportDescriptors = reportDescriptors;
            this.transformerMap = new HashMap<>();
            if (transformers != null){
                for (GATKReportTableTransformer t : transformers){
                    this.transformerMap.put(t.getEvalModuleName(), t);
                }
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

        public void configureWalker(VariantQC variantQC) {
            this.walker = new VariantEvalChild(variantQC, this);
        }
    }

    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        if (pedigreeFile != null)
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));

        Collection<String> samples = getHeaderForVariants().getSampleNamesInOrder();
        if (samples != null) {
            sampleDBBuilder.addSamplesFromSampleNames(samples);
        }

        return sampleDBBuilder.getFinalSampleDB();
    }

    //TODO: add getter to VariantEvalUtils
    private Map<String, Class<? extends VariantStratifier>> getEvaluationClassMap() {
        Map<String, Class<? extends VariantStratifier>> stratifierClasses = new HashMap<>();
        Reflections reflectionsStrat = new Reflections("org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications");
        Set<Class<? extends VariantStratifier>> allClasses = reflectionsStrat.getSubTypesOf(VariantStratifier.class);
        for (Class<? extends VariantStratifier> clazz : allClasses) {
            stratifierClasses.put(clazz.getSimpleName(), clazz);
        }

        return stratifierClasses;
    }
}
