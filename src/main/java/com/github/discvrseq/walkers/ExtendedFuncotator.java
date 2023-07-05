package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvValidationException;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This is an extension of the GATK Funcotator tool. The primary purpose of this tool is to produce a VCF output where specific fields from Funcotator are output
 * as discrete INFO field annotations, rather than Funcotator's default, which involves a single INFO field where the annotations are concatenated together.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar ExtendedFuncotator \
 *     -V input.vcf.gz \
 *     -R myFasta.fasta \
 *     --ref-version hg19 \
 *     --data-sources-path ./path/to/dataSource \
 *     -cf configFile.txt \
 *     -O output.annotated.vcf.gz
 * </pre>
 *
 * The config file must specify the list of Funcotator fields to include, and the VCF info attributes for each. The fields ID, Number, Type, and Description must match valid VCF header values.
 * <pre>
 *    DataSource	SourceField	ID	Number	Type	Description
 *    testSource	OtherField	f1	UNBOUNDED	String	This is field1
 *    testSource	ExpectedResult	f2	UNBOUNDED	String	This is field2
 *    vcfTestSource	AA	f3	A	Float	This is field3
 *    vcfTestSource	BB	f4	UNBOUNDED	String	This is field4
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Create functional annotations on given variants cross-referenced by a given set of data sources. Based on GATK Funcotator, but with expanded output options.\n",
        oneLineSummary = "Functional Annotator",
        programGroup = DiscvrSeqProgramGroup.class
)
@DocumentedFeature
public class ExtendedFuncotator extends Funcotator {
    private static final Logger logger = LogManager.getLogger(ExtendedFuncotator.class);

    @Argument(doc = "A TSV file specifying the Funcotator fields to extract, and their VCF INFO field annotation information.", fullName = "config-file", shortName = "cf", optional = false)
    public GATKPath configFile = null;

    @Argument(doc = "If specified, a list of all fields present in the data sources but not included in the output will be printed, and the tool will immediately exit.", fullName = "print-missing-fields", shortName = "pmf", optional = true)
    public boolean printMissingFields = false;

    @Override
    public void onTraversalStart() {
        // This will be ignored anyway
        getArguments().outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;

        addOutputVCFCommandLine = false;

        // Get our overrides for annotations:
        final LinkedHashMap<String, String> annotationDefaultsMap = FuncotatorEngine.splitAnnotationArgsIntoMap(getArguments().annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = FuncotatorEngine.splitAnnotationArgsIntoMap(getArguments().annotationOverrides);

        // Get the header for our variants:
        final VCFHeader vcfHeader = getHeaderForVariants();

        final Set<String> finalUserTranscriptIdSet = FuncotatorEngine.processTranscriptList(getArguments().userTranscriptIdSet);
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths( getArguments().referenceVersion,  getArguments().dataSourceDirectories);

        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                annotationOverridesMap,
                getArguments().transcriptSelectionMode,
                finalUserTranscriptIdSet,
                this,
                getArguments().lookaheadFeatureCachingInBp,
                new FlankSettings(getArguments().fivePrimeFlankSize, getArguments().threePrimeFlankSize),
                false,
                getArguments().minNumBasesForValidSegment
        );


        // Read config, etc:
        List<VcfHeaderDescriptor> fields = new ArrayList<>();
        IOUtil.assertFileIsReadable(configFile.toPath());
        try (CSVReader reader = CsvUtils.getTsvReader(configFile.toPath().toFile()))
        {
            String[] line;
            List<String> header = new ArrayList<>();
            int idx = 0;
            while ((line = reader.readNext()) != null)
            {
                idx++;

                if (idx == 1)
                {
                    header.addAll(Arrays.stream(line).map(String::toLowerCase).toList());
                    continue;
                }

                fields.add(new VcfHeaderDescriptor(line, header));
            }
        }
        catch (IOException | CsvValidationException e)
        {
            throw new GATKException("Unable to read config file", e);
        }

        fields.forEach(x -> vcfHeader.addMetaDataLine(x.toInfoLine()));
        List<VCFInfoHeaderLine> headerLines = new ArrayList<>(vcfHeader.getInfoHeaderLines());

        funcotatorEngine = new FuncotatorEngine(
                getArguments(),
                getSequenceDictionaryForDrivingVariants(),
                VcfFuncotationMetadata.create(headerLines),
                dataSourceFuncotationFactories
        );

        if (printMissingFields) {
            Map<String, List<String>> allFields = new HashMap<>();

            dataSourceFuncotationFactories.stream().forEach(f -> {
                allFields.put(f.getName(), new ArrayList<>(f.getSupportedFuncotationFields()));
            });

            fields.forEach(f -> {
                if (!allFields.containsKey(f.dataSource)) {
                    logger.warn("Requested field from non-existing data source: " + f.dataSource);
                    return;
                }

                String target = f.dataSource + "_" + f.sourceField;
                allFields.get(f.dataSource).remove(target);
            });

            allFields.keySet().forEach(f -> {
                if (allFields.get(f).isEmpty()) {
                    return;
                }

                logger.info("The following fields are present but not used from: " + f);
                allFields.get(f).forEach(logger::info);
            });

            System.exit(0);
        }

        // Create our output renderer:
        logger.info("Creating output: " + getArguments().outputFile.toURI());
        outputRenderer = new ExtendedVcfOutputRenderer(
                this.createVCFWriter(getArguments().outputFile),
                funcotatorEngine.getFuncotationFactories(),
                vcfHeader,
                annotationDefaultsMap,
                annotationOverridesMap,
                getDefaultToolVCFHeaderLines(),
                getArguments().excludedFields,
                this.getVersion(),
                fields
        );
    }

    private static class VcfHeaderDescriptor
    {
        final String dataSource;
        final String sourceField;
        final String id;
        final VCFHeaderLineCount count;
        final VCFHeaderLineType type;
        final String description;

        public VcfHeaderDescriptor(String[] line, List<String> headerFields)
        {
            dataSource = getField("DataSource", line, headerFields);
            sourceField = getField("SourceField", line, headerFields);
            id = getField("ID", line, headerFields);

            String numberString = getField("Number", line, headerFields);
            try {
                count = VCFHeaderLineCount.valueOf(numberString);
            }
            catch (IllegalArgumentException e) {
                throw new GATKException("Unknown value for Number columns: " + numberString);
            }

            String typeString = getField("Type", line, headerFields);
            try {
                type = VCFHeaderLineType.valueOf(typeString);
            }
            catch (IllegalArgumentException e)
            {
                throw new GATKException("Unknown value for type: " + typeString);
            }

            description = getField("Description", line, headerFields);
        }

        private String getField(String fieldName, String[] line, List<String> headerFields)
        {
            if (!headerFields.contains(fieldName.toLowerCase()))
            {
                throw new GATKException("Config file should contain the column: " + fieldName);
            }

            int idx = headerFields.indexOf(fieldName.toLowerCase());
            if (idx > line.length)
            {
                throw new GATKException("Missing value for column: " + fieldName);
            }

            return line[idx];
        }

        public VCFInfoHeaderLine toInfoLine()
        {
            return new VCFInfoHeaderLine(id, count, type, description);
        }
    }

    protected void enqueueAndHandleVariant(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForVariant(variant, referenceContext, featureContext);

        // At this point there is only one transcript ID in the funcotation map if canonical or best effect are selected
        outputRenderer.write(variant, funcotationMap);
    }

    private static final class ExtendedVcfOutputRenderer extends VcfOutputRenderer {
        private final VariantContextWriter vcfWriter;

        private final List<VcfHeaderDescriptor> fieldsToOutput;

        private final Set<String> keysEncountered = new HashSet<>();

        public ExtendedVcfOutputRenderer(final VariantContextWriter vcfWriter,
                                 final List<DataSourceFuncotationFactory> dataSources,
                                 final VCFHeader existingHeader,
                                 final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                 final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                 final Set<VCFHeaderLine> defaultToolVcfHeaderLines,
                                 final Set<String> excludedOutputFields,
                                 final String toolVersion,
                                 final List<VcfHeaderDescriptor> fieldsToOutput) {
            super(vcfWriter, dataSources, existingHeader, unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, defaultToolVcfHeaderLines, excludedOutputFields, toolVersion);
            this.vcfWriter = vcfWriter;

            this.fieldsToOutput = fieldsToOutput;
        }

        @Override
        public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

            // Create a new variant context builder:
            final VariantContextBuilder variantContextOutputBuilder = new VariantContextBuilder(variant);

            // Add our new annotations:
            for (VcfHeaderDescriptor vd : fieldsToOutput) {
                final String target = vd.dataSource + "_" + vd.sourceField;
                final Map<Allele, List<String>> valMap = new HashMap<>();

                for (final String txId : txToFuncotationMap.getTranscriptList()) {
                    final List<Funcotation> funcotations = txToFuncotationMap.get(txId);
                    List<Funcotation> funcotationList = funcotations.stream()
                            .filter(f -> f.getFieldNames().size() > 0)
                            .filter(f -> f.getDataSourceName().equalsIgnoreCase(vd.dataSource))
                            .filter(f -> f.getFieldNames().contains(target))
                            .toList();

                    if (funcotationList.isEmpty()) {
                        continue;
                    }

                    funcotationList.forEach(f -> {
                        if (!valMap.containsKey(f.getAltAllele())) {
                            valMap.put(f.getAltAllele(), new ArrayList<>());
                        }

                        if (f.getField(target) != null && !f.getField(target).isEmpty()) {
                            valMap.get(f.getAltAllele()).add(f.getField(target));
                        }
                    });
                }

                List<List<String>> toAdd;
                switch (vd.count) {
                    case A -> {
                        toAdd = variant.getAlternateAlleles().stream().map(a -> valMap.getOrDefault(a, Collections.emptyList())).toList();
                        if (toAdd.stream().anyMatch(x -> x.size() > 1)) {
                            throw new GATKException("Expected a 0-1 values per allele for annotation: " + vd.sourceField + ". Problem at site: " + variant.toStringWithoutGenotypes());
                        }

                        // We expect only one value per allele, so grab the first:
                        if (!toAdd.isEmpty() && !toAdd.stream().map(x -> x.isEmpty() ? "" : x.get(0)).allMatch(String::isEmpty)) {
                            String val = toAdd.stream().map(x -> x.isEmpty() ? "" : x.get(0)).collect(Collectors.joining(","));
                            if (!val.isEmpty()) {
                                variantContextOutputBuilder.attribute(vd.id, val);
                            }
                        }
                    }
                    case R -> {
                        toAdd = variant.getAlleles().stream().map(a -> valMap.getOrDefault(a, Collections.emptyList())).toList();
                        if (toAdd.stream().anyMatch(x -> x.size() > 1)) {
                            throw new GATKException("Expected a 0-1 values per allele for annotation: " + vd.sourceField + ". Problem at site: " + variant.toStringWithoutGenotypes());
                        }

                        // We expect only one value per allele, so grab the first:
                        if (!toAdd.isEmpty() && !toAdd.stream().map(x -> x.isEmpty() ? "" : x.get(0)).allMatch(String::isEmpty)) {
                            String val = toAdd.stream().map(x -> x.isEmpty() ? "" : x.get(0)).collect(Collectors.joining(","));
                            if (!val.isEmpty()) {
                                variantContextOutputBuilder.attribute(vd.id, val);
                            }
                        }
                    }
                    case INTEGER, UNBOUNDED -> {
                        // The assumption is that if the annotation applies to every ALT allele, it will get annotated on each
                        // We therefore take the unique values annotated to each allele, which should produce one list
                        toAdd = valMap.values().stream().distinct().toList();
                        if (toAdd.size() > 1 && vd.count != VCFHeaderLineCount.UNBOUNDED) {
                            String values = valMap.keySet().stream().map(key -> key + ":" + StringUtils.join(valMap.get(key), ",")).collect(Collectors.joining("; "));
                            throw new GATKException("Expected a one set of values per site for annotation: " + vd.dataSource + " / " + vd.sourceField + ". Problem at site: " + variant.toStringWithoutGenotypes() + ", with annotation: " + values);
                        }

                        if (!toAdd.isEmpty() && !toAdd.get(0).isEmpty()) {
                            String val = toAdd.get(0).stream().filter(x -> !x.isEmpty()).collect(Collectors.joining(","));
                            if (!val.isEmpty()) {
                                variantContextOutputBuilder.attribute(vd.id, val);
                            }
                        }
                    }
                }
            }

            // Add the genotypes from the variant:
            variantContextOutputBuilder.genotypes(variant.getGenotypes());

            // Render and add our VCF line:
            vcfWriter.add(variantContextOutputBuilder.make());
        }
    }
}
