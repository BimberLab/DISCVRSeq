package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import java.io.File;
import java.util.*;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
        summary = "Annotate a VCF with clinical variants using ClinVar vcf_2.0",
        oneLineSummary = "Annotate a VCF with ClinVar vcf_2.0",
        programGroup = DiscvrSeqProgramGroup.class
)
public class ClinvarAnnotator extends VariantWalker {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    @Argument(doc="Clinvar VCF", fullName = "clinvar", shortName = "clinvar", optional = false)
    public List<FeatureInput<VariantContext>> clinvarVariants = null;

    private VariantContextWriter writer = null;

    private List<VCFInfoHeaderLine> HEADER_LINES = Arrays.asList(
            new VCFInfoHeaderLine("CLN_ALLELE", VCFHeaderLineCount.R, VCFHeaderLineType.Character, "Alternate alleles from Clinvar"),
            new VCFInfoHeaderLine("CLN_ALLELEID", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "the ClinVar Allele ID"),
            new VCFInfoHeaderLine("CLN_DN", VCFHeaderLineCount.R, VCFHeaderLineType.String, "ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB"),
            new VCFInfoHeaderLine("CLN_DNINCL", VCFHeaderLineCount.R, VCFHeaderLineType.String, "For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB"),
            new VCFInfoHeaderLine("CLN_DISDB", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN"),
            new VCFInfoHeaderLine("CLN_DISDBINCL", VCFHeaderLineCount.R, VCFHeaderLineType.String,"For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN"),
            new VCFInfoHeaderLine("CLN_HGVS", VCFHeaderLineCount.R, VCFHeaderLineType.String,"For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN"),
            new VCFInfoHeaderLine("CLN_REVSTAT", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Top-level (primary assembly, alt, or patch) HGVS expression."),
            new VCFInfoHeaderLine("CLN_SIG", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Clinical significance for this single variant"),
            new VCFInfoHeaderLine("CLN_SIGINCL", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance."),
            new VCFInfoHeaderLine("CLN_VC", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Variant type"),
            new VCFInfoHeaderLine("CLN_VCSO", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Sequence Ontology id for variant type"),
            new VCFInfoHeaderLine("CLN_VI", VCFHeaderLineCount.R, VCFHeaderLineType.String,"the variant's clinical sources reported as tag-value pairs of database and variant identifier"),
            new VCFInfoHeaderLine("CLN_DBVARID", VCFHeaderLineCount.R, VCFHeaderLineType.String, "nsv accessions from dbVar for the variant"),
            new VCFInfoHeaderLine("CLN_GENEINFO", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)"),
            new VCFInfoHeaderLine("CLN_MC", VCFHeaderLineCount.R, VCFHeaderLineType.String,"comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence"),
            new VCFInfoHeaderLine("CLN_ORIGIN", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other"),
            new VCFInfoHeaderLine("CLN_RS", VCFHeaderLineCount.R, VCFHeaderLineType.String,"dbSNP ID (i.e. rs number)"),
            new VCFInfoHeaderLine("CLN_SSR", VCFHeaderLineCount.R, VCFHeaderLineType.String,"Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other")
    );

    @Override
    public void onTraversalStart() {
        VCFHeader header = new VCFHeader(getHeaderForVariants());
        for (VCFInfoHeaderLine line : HEADER_LINES){
            header.addMetaDataLine(line);
        }

        writer = createVCFWriter(out);
        writer.writeHeader(header);
    }

    int max = 0;
    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        max = variant.getStart();
        Map<Allele, Map<String, String>> annotationMap = new HashMap<>();
        List<VariantContext> matches = featureContext.getValues(clinvarVariants);
        if (matches.isEmpty()){
            return;
        }

        boolean foundHit = false;
        for (VariantContext vc : matches){
            if (vc.getAlternateAlleles().size() != 1){
                throw new IllegalStateException("More than 1 alternate allele found, please use new clinvar vcf_2.0 with 1 alt allele per position");
            }

            if (!variant.getReference().equals(vc.getReference())) {
                //throw new IllegalStateException("different references");
                logger.warn("Reference not equal(chr:startPos:endPos Ref[ALTs]) to ClinVar's:\t" + variant.getContig() + ":" + variant.getStart() + ":" + variant.getEnd() + "\t" + variant.getReference() + variant.getAlternateAlleles() + "\t" + vc.getContig() + ":" + vc.getStart() + ":" + vc.getEnd() + "\t" + vc.getReference() + vc.getAlternateAlleles());
                continue;
            }

            if (variant.getStart() != vc.getStart()) {
                //throw new IllegalStateException("different start positions");
                logger.warn("Start position not equal(chr:startPos:endPos Ref[ALTs]) to ClinVar's:\t" + variant.getContig() + ":" + variant.getStart() + ":" + variant.getEnd() + "\t" + variant.getReference() + variant.getAlternateAlleles() + "\t" + vc.getContig() + ":" + vc.getStart() + ":" + vc.getEnd() + "\t" + vc.getReference() + vc.getAlternateAlleles());
                continue;
            }

            if (variant.getEnd() != vc.getEnd()) {
                logger.warn("End position not equal(chr:startPos:endPos Ref[ALTs]) to ClinVar's:\t" + variant.getContig() + ":" + variant.getStart() + ":" + variant.getEnd() + "\t" + variant.getReference() + variant.getAlternateAlleles() + "\t" + vc.getContig() + ":" + vc.getStart() + ":" + vc.getEnd() + "\t" + vc.getReference() + vc.getAlternateAlleles());
                continue;
            }

            Allele cvAllele = vc.getAlternateAllele(0);
            for (Allele alt : variant.getAlternateAlleles()){
                if (cvAllele.equals(alt)){
                    //gather annotations, add to map by allele
                    foundHit = true;
                    annotationMap.put(alt, transferAnnotations(vc, alt, vcb));
                }
            }
        }

        if (!foundHit){
            return;
        }

        for (VCFInfoHeaderLine line : HEADER_LINES){
            List<String> sb = new ArrayList<>();

            int nonNull = 0;
            for (Allele a : vcb.getAlleles()){
                if (annotationMap.containsKey(a) && annotationMap.get(a).get(line.getID()) != null ) {
                    sb.add(annotationMap.get(a).get(line.getID()));
                    nonNull++;
                }
                else {
                    sb.add("");
                }
            }

            //Only include annotations with values
            if (nonNull > 0){
                if (!sb.isEmpty()){
                    vcb.attribute(line.getID(), sb.stream().collect(Collectors.joining(",")));
                }
            }
        }

        writer.add(vcb.make());
    }

    /**
     * Transfer Annotations from ClinVar VCF to new VCF
     * @param source    ClinVar VCF
     * @param alt   alternate allele from ClinVar
     * @param vcb   new VCF to annotate
     * @return
     */
    private Map<String, String> transferAnnotations (VariantContext source, Allele alt, VariantContextBuilder vcb){
        //TODO: extract each annotation from source (for that specific alt allele)
        //determine index of the alt allele in the target variant
        //transfer that annotation to appropriate position in the target

        //examples:
        //Ref: A, Alt: T,C
        //Clinvar: C,G
        //Clinvar annotations are comma separated arrays, with each position matching the corresponding ALT allele
        //The allele C matches, so transfer matching clinvar annotation to target VCF.

        Map<String, String> annotations = new HashMap<>();

        annotations.put("CLN_ALLELE", alt.toString());
        annotations.put("CLN_ALLELEID", annotateValue(source,"ALLELEID"));
        annotations.put("CLN_DN", annotateValue(source,"CLNDN"));
        annotations.put("CLN_DNINCL", annotateValue(source, "CLNDNINCL"));
        annotations.put("CLN_DISDB", annotateValue(source, "CLNDISDB"));
        annotations.put("CLN_DISDBINCL", annotateValue(source,"CLN_DISDBINCL"));
        annotations.put("CLN_HGVS", annotateValue(source,"CLN_HGVS"));
        annotations.put("CLN_REVSTAT", annotateValue(source,"CLN_REVSTAT"));
        annotations.put("CLN_SIG", annotateValue(source,"CLN_SIG"));
        annotations.put("CLN_SIGINCL", annotateValue(source,"CLN_SIGINCL"));
        annotations.put("CLN_VC", annotateValue(source,"CLN_VC"));
        annotations.put("CLN_VCSO", annotateValue(source,"CLN_VCSO"));
        annotations.put("CLN_VI", annotateValue(source,"CLN_VI"));
        annotations.put("CLN_GENEINFO", annotateValue(source,"CLN_GENEINFO"));
        annotations.put("CLN_MC", annotateValue(source,"CLN_MC"));
        annotations.put("CLN_ORIGIN", annotateValue(source,"CLN_ORIGIN"));
        annotations.put("CLN_RS", annotateValue(source,"CLN_RS"));
        annotations.put("CLN_SSR", annotateValue(source,"CLN_SSR"));

        return annotations;
    }

    /**
     * Obtain value for given ClinVar VCF annotation ID
     * @param source    ClinVar VCF
     * @param ID   ClinVar VCF Annotation ID
     * @return
     */
    private String annotateValue(VariantContext source, String ID){
        if (source.getAttribute(ID) == null){
            return null;
        }

        //Join array annotations by | character
        return source.getAttribute(ID) instanceof Collection ? StringUtils.join(source.getAttributeAsList(ID), "|") : source.getAttributeAsString(ID, "");
    }

    @Override
    public Object onTraversalSuccess() {
        logger.warn("max: " + max);
        return super.onTraversalSuccess();
    }

    /**
     * Close out the new variants file.
     */
    @Override
    public void closeTool() {
        if (writer != null) {
            writer.close();
        }
    }
}
