package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

import java.io.File;
import java.util.*;

/**
 *
 * This tool can be used to remove specific INFO or genotype annotations from a VCF (while retaining others), or to completely clear annotations.
 *
 * <h3>Usage examples:</h3>
 * <h3>Discard all annotations except AF:</h3>
 * <pre>
 *   java -jar DISCVRSeq.jar RemoveAnnotations \
 *     -variant input.vcf \
 *     -A AF \
 *     -O output.vcf
 * </pre>
 *
 * <h3>Discard only the MAF annotation:</h3>
 * <pre>
 *   java -jar DISCVRSeq.jar RemoveAnnotations \
 *     -variant input.vcf \
 *     -XA MAF \
 *     -O output.vcf
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Removes specific site or genotype annotations from a VCF",
        oneLineSummary = "Removes specific site or genotype annotations from a VCF",
        programGroup = DiscvrSeqProgramGroup.class
)
public class RemoveAnnotations extends VariantWalker {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    @Argument(fullName="annotationToKeep", shortName="A", doc="List the specific INFO field annotations to retain (by their keys, such as AC, AF, etc).  If specified, all other INFO field annotations will be removed.", optional=true)
    protected List<String> annotationToKeep = new ArrayList<>();

    @Argument(fullName="annotationToRemove", shortName="XA", doc="One or more specific INFO field annotations to remove (by their keys, such as AC, AF, etc).", optional=true)
    protected List<String> annotationsToExclude = new ArrayList<>();

    @Argument(fullName="genotypeAnnotationToKeep", shortName="GA", doc="List the specific genotype (format field) annotations to retain", optional=true)
    protected List<String> genotypeAnnotationToKeep = new ArrayList<>();

    @Argument(fullName="genotypeAnnotationToRemove", shortName="XGA", doc="One or more specific genotype (format field) annotations to remove.", optional=true)
    protected List<String> genotypeAnnotationsToExclude = new ArrayList<>();

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered sites", optional=true)
    protected boolean excludeFiltered = true;

    @Argument(fullName="retainExtraHeaderLines", shortName="rh", doc="If provided, additional header lines (metadata, etc) will be retained.", optional=true)
    protected boolean retainExtraHeaderLines = false;

    @Argument(fullName="clearGenotypeFilter", shortName="cgf", doc="Clear the filter field on all genotypes.  This executes after setFilteredGTToNoCall", optional=true)
    protected boolean clearGTfilter = true;

    @Argument(fullName="setFilteredGTToNoCall", shortName="sgf", doc="Sets filtered genotypes to no-call.  ", optional=true)
    protected boolean setFilteredGTToNoCall = true;

    @Argument(fullName="sitesOnly", shortName="sitesOnly", doc="Omit samples and genotypes from the output VCF.  ", optional=true)
    protected boolean sitesOnly = false;

    private VCFHeader header;
    private Set<String> allowableInfoKeys;
    private Set<String> allowableFormatKeys;

    private List<String> formatFields = Arrays.asList(
            VCFConstants.GENOTYPE_KEY,
            VCFConstants.DEPTH_KEY,
            VCFConstants.GENOTYPE_PL_KEY,
            VCFConstants.GENOTYPE_QUALITY_KEY,
            VCFConstants.GENOTYPE_ALLELE_DEPTHS
    );

    private VariantContextWriter writer = null;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        VCFHeader initialHeader = new VCFHeader(getHeaderForVariants());

        //retain only specified header lines
        Set<VCFHeaderLine> headerLines = new HashSet<>();
        int skippedHeaderLines = 0;

        //strip info annotations
        for (VCFInfoHeaderLine line : initialHeader.getInfoHeaderLines()) {
            skippedHeaderLines += inspectAnnotation(headerLines, line, annotationToKeep, annotationsToExclude);
        }

        //strip format fields
        for (VCFFormatHeaderLine line : initialHeader.getFormatHeaderLines()) {
            //special-case standard genotype field
            if (!sitesOnly && formatFields.contains(line.getID())){
                headerLines.add(line);
                continue;
            }

            skippedHeaderLines += inspectAnnotation(headerLines, line, genotypeAnnotationToKeep, genotypeAnnotationsToExclude);
        }

        //strip filters, if selected
        if (!excludeFiltered){
            headerLines.addAll(initialHeader.getFilterLines());
        }

        //metadata
        if (retainExtraHeaderLines){
            headerLines.addAll(initialHeader.getOtherHeaderLines());
        }
        else {
            skippedHeaderLines += initialHeader.getOtherHeaderLines().size();
        }

        headerLines.addAll(initialHeader.getContigLines());

        logger.info("total header lines skipped: " + skippedHeaderLines);

        header = new VCFHeader(headerLines, (sitesOnly ? Collections.emptyList() : initialHeader.getSampleNamesInOrder()));

        allowableInfoKeys = new HashSet<>();
        //NOTE: if lambdas are used here, the walker will not be picked up by PluginManager
        // http://gatkforums.broadinstitute.org/gatk/discussion/comment/38892#Comment_38892
        for (VCFInfoHeaderLine line : header.getInfoHeaderLines()){
            allowableInfoKeys.add(line.getID());
        }

        allowableFormatKeys = new HashSet<>();
        for (VCFFormatHeaderLine line : header.getFormatHeaderLines()){
            allowableFormatKeys.add(line.getID());
        }

        writer = createVCFWriter(out);
        writer.writeHeader(header);
    }

    private int inspectAnnotation(Set<VCFHeaderLine> headerLines, VCFCompoundHeaderLine line, List<String> annotationToKeep, List<String> annotationsToExclude){
        if (!annotationToKeep.isEmpty() && !annotationToKeep.contains(line.getID())){
            return 1;
        }

        if (annotationsToExclude.contains(line.getID())){
            return 1;
        }

        headerLines.add(line);
        return 0;
    }

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (excludeFiltered && vc.isFiltered()){
            return;
        }

        VariantContextBuilder vcb = new VariantContextBuilder(vc);

        //retain only allowable info fields
        Set<String> keys = new HashSet<>(vc.getAttributes().keySet());
        keys.removeAll(allowableInfoKeys);
        vcb.rmAttributes(new ArrayList<>(keys));

        //now genotypes:
        if (sitesOnly){
            vcb.noGenotypes();
        }
        else {
            GenotypesContext ctx = GenotypesContext.copy(vc.getGenotypes());
            Set<String> sampleNames = new HashSet<>(ctx.getSampleNames());
            for (String sn : sampleNames){
                Genotype g = ctx.get(sn);
                GenotypeBuilder gb = new GenotypeBuilder(g);

                if (setFilteredGTToNoCall && g.isFiltered()){
                    gb.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
                    gb.noAD();
                    gb.noDP();
                    gb.noGQ();
                    gb.noPL();
                    gb.noAttributes();
                }

                if (clearGTfilter){
                    gb.unfiltered();
                }

                //retain only allowable info fields
                Map<String, Object> gtAttributes = new HashMap<>(g.getExtendedAttributes());
                gtAttributes.keySet().retainAll(allowableFormatKeys);
                gb.noAttributes();
                gb.attributes(gtAttributes);

                ctx.replace(gb.make());
            }
            vcb.genotypes(ctx);
        }

        writer.add(vcb.make());
    }

    /**
     * Closes out the new variants file.
     */
    @Override
    public void closeTool() {
        if (writer != null) {
            writer.close();
        }
    }
}
