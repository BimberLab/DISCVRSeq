package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.*;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Create functional annotations on given variants cross-referenced by a given set of data sources.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar ExtendedFuncotator \
 *     -V input.vcf.gz \
 *     -O output.annotated.vcf.gz
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Create functional annotations on given variants cross-referenced by a given set of data sources.\n" +
                "A GATK functional annotation tool (similar functionality to Oncotator).",
        oneLineSummary = "Functional Annotator",
        programGroup = DiscvrSeqProgramGroup.class
)
@DocumentedFeature
public class ExtendedFuncotator extends Funcotator {
    private static final Logger logger = LogManager.getLogger(ExtendedFuncotator.class);

    @Override
    public void onTraversalStart() {

        // Get our overrides for annotations:
        final LinkedHashMap<String, String> annotationDefaultsMap = FuncotatorEngine.splitAnnotationArgsIntoMap(getArguments().annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = FuncotatorEngine.splitAnnotationArgsIntoMap(getArguments().annotationOverrides);

        // Get the header for our variants:
        final VCFHeader vcfHeader = getHeaderForVariants();

        // Create our output renderer:
        logger.info("Creating a " + getArguments().outputFormatType + " file for output: " + getArguments().outputFile.toURI());
        outputRenderer = new ExtendedVcfOutputRenderer(
                this.createVCFWriter(getArguments().outputFile),
                funcotatorEngine.getFuncotationFactories(),
                vcfHeader,
                annotationDefaultsMap,
                annotationOverridesMap,
                getDefaultToolVCFHeaderLines(),
                getArguments().excludedFields,
                this.getVersion()
        );
    }

    protected void enqueueAndHandleVariant(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForVariant(variant, referenceContext, featureContext);

        // At this point there is only one transcript ID in the funcotation map if canonical or best effect are selected
        outputRenderer.write(variant, funcotationMap);
    }

    private static final class ExtendedVcfOutputRenderer extends VcfOutputRenderer {
        private final VariantContextWriter vcfWriter;

        public ExtendedVcfOutputRenderer(final VariantContextWriter vcfWriter,
                                 final List<DataSourceFuncotationFactory> dataSources,
                                 final VCFHeader existingHeader,
                                 final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                                 final LinkedHashMap<String, String> unaccountedForOverrideAnnotations,
                                 final Set<VCFHeaderLine> defaultToolVcfHeaderLines,
                                 final Set<String> excludedOutputFields,
                                 final String toolVersion) {
            super(vcfWriter, dataSources, existingHeader, unaccountedForDefaultAnnotations, unaccountedForOverrideAnnotations, defaultToolVcfHeaderLines, excludedOutputFields, toolVersion);
            this.vcfWriter = vcfWriter;
        }

        @Override
        public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {

            // Create a new variant context builder:
            final VariantContextBuilder variantContextOutputBuilder = new VariantContextBuilder(variant);

            final StringBuilder funcotatorAnnotationStringBuilder = new StringBuilder();

            // Get the old VCF Annotation field and append the new information to it:
            final Object existingAnnotation = variant.getAttribute(FUNCOTATOR_VCF_FIELD_NAME, null);
            final List<String> existingAlleleAnnotations;
            if ( existingAnnotation != null) {
                existingAlleleAnnotations = Utils.split(existingAnnotation.toString(), ',');
            }
            else {
                existingAlleleAnnotations = Collections.emptyList();
            }

            // Go through each allele and add it to the writer separately:
            final List<Allele> alternateAlleles = variant.getAlternateAlleles();
            for ( int alleleIndex = 0; alleleIndex < alternateAlleles.size() ; ++alleleIndex ) {

                final Allele altAllele = alternateAlleles.get(alleleIndex);

                if ( alleleIndex < existingAlleleAnnotations.size() ) {
                    funcotatorAnnotationStringBuilder.append( existingAlleleAnnotations.get(alleleIndex) );
                    funcotatorAnnotationStringBuilder.append(FIELD_DELIMITER);
                }

                for (final String txId : txToFuncotationMap.getTranscriptList()) {
                    funcotatorAnnotationStringBuilder.append(START_TRANSCRIPT_DELIMITER);
                    final List<Funcotation> funcotations = txToFuncotationMap.get(txId);
                    final Funcotation manualAnnotationFuncotation = createManualAnnotationFuncotation(altAllele);

                    funcotatorAnnotationStringBuilder.append(
                            Stream.concat(funcotations.stream(), Stream.of(manualAnnotationFuncotation))
                                    .filter(f -> f.getAltAllele().equals(altAllele))
                                    .filter(f -> f.getFieldNames().size() > 0)
                                    .filter(f -> !f.getDataSourceName().equals(FuncotatorConstants.DATASOURCE_NAME_FOR_INPUT_VCFS))
                                    .map(VcfOutputRenderer::adjustIndelAlleleInformation)
                                    .map(f -> FuncotatorUtils.renderSanitizedFuncotationForVcf(f, finalFuncotationFieldNames))
                                    .collect(Collectors.joining(FIELD_DELIMITER))
                    );

                    funcotatorAnnotationStringBuilder.append(END_TRANSCRIPT_DELIMITER + ALL_TRANSCRIPT_DELIMITER);
                }
                // We have a trailing "#" - we need to remove it:
                funcotatorAnnotationStringBuilder.deleteCharAt(funcotatorAnnotationStringBuilder.length()-1);
                funcotatorAnnotationStringBuilder.append(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
            }

            // We have a trailing "," - we need to remove it:
            funcotatorAnnotationStringBuilder.deleteCharAt(funcotatorAnnotationStringBuilder.length()-1);

            // Add our new annotation:
            variantContextOutputBuilder.attribute(FUNCOTATOR_VCF_FIELD_NAME, funcotatorAnnotationStringBuilder.toString());

            // Add the genotypes from the variant:
            variantContextOutputBuilder.genotypes( variant.getGenotypes() );

            // Render and add our VCF line:
            vcfWriter.add( variantContextOutputBuilder.make() );
        }
    }
}
