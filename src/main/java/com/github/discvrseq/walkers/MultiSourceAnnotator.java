package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import com.github.discvrseq.walkers.annotator.Impact;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;

/**
 * This is a fairly specialized tool, designed to take the VCFs annotated with ClinvarAnnotator (from DISCVR-seq Toolkit) and <a href="https://www.hgsc.bcm.edu/software/cassandra">Cassandra</a>,
 * and transfer select annotations from those VCFs to the input VCF.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar MultiSourceAnnotator \
 *     -cv clinvar.vcf.gz \
 *     -lr liftoverFailures.vcf \
 *     -c cassandraVcf.vcf.gz \
 *     -O output.vcf.gz
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This is a fairly specialized tool, designed to take the VCFs annotated with ClinvarAnnotator and Cassandra, and transfer select annotations from those VCFs to the input VCF.",
        oneLineSummary = "Transfers annotations from ClinvarAnnotator and Cassandra VCFs to a source VCF",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class MultiSourceAnnotator extends VariantWalker {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(doc="Clinvar Annotated VCF", fullName = "clinvar", shortName = "cv", optional = true)
    public FeatureInput<VariantContext> clinvarVariants = null;

    @Argument(doc="Liftover Reject VCF", fullName = "liftoverReject", shortName = "lr", optional = true)
    public FeatureInput<VariantContext> liftoverRejectVariants = null;

    @Argument(doc="Cassandra Annotated VCF", fullName = "cassandra", shortName = "c", optional = true)
    public FeatureInput<VariantContext> cassandraVariants = null;

    @Argument(doc="SnpSift Annotated VCF", fullName = "snpsift", shortName = "ss", optional = true)
    public FeatureInput<VariantContext> snpSiftVariants = null;

    @Argument(doc="Funcotator Annotated VCF", fullName = "funcotator", shortName = "f", optional = true)
    public FeatureInput<VariantContext> funcotatorVariants = null;

    private VariantContextWriter writer;

    private final List<String> CLINVAR_INFO = Arrays.asList(
            "CLN_ALLELE",
            "CLN_ALLELEID",
            "CLN_DBVARID",
            "CLN_DISDB",
            "CLN_DISDBINCL",
            "CLN_DN",
            "CLN_DNINCL",
            "CLN_GENEINFO",
            "CLN_HGVS",
            "CLN_MC",
            "CLN_ORIGIN",
            "CLN_REVSTAT",
            "CLN_RS",
            "CLN_SIG",
            "CLN_SIGINCL",
            "CLN_SSR",
            "CLN_VC",
            "CLN_VCSO",
            "CLN_VI"
    );

    private final List<String> FUNCOTATOR_INFO = Arrays.asList(
            "ACMG_Disease",
            "ACMGLMM_LOF",
            "ACMGLMM_MOI",
            "ACMGLMM_GN",
            "CADD_Score",
            "CGC_Synd",
            "CGC_GeneID",
            "CGC_Mut _Type",
            "CGC_GeneID",
            "CGC_Tissue",
            "CGC_Germline",
            "CGC_Somatic",
            "CAPICE",
            "ClinPredScore",
            "CLN_ALLELEID",
            "CLN_DN",
            "CLN_Name",
            "CLNVarReviewStatus",
            "CLN_CS",
            "CLN_CSC",
            "DBSCSNV",
            "AFby1000KG",
            "CommonSNP",
            "dbSNPBUildID",
            "RS",
            "DNASECLUST",
            "ERB_TYPE",
            "Fam_Synd",
            "FANTOM_ENHNCRID",
            "FANTOM_SCORE",
            "FANTOM_STRND",
            "FANTOM_TFBS_SCORE",
            "FANTOM_TFBS_STRAND",
            "FANTOM_TF",
            "F_C",
            "F_CG",
            "F_NC",
            "F_NCG",
            "FS2",
            "AF_gnomAD",
            "gnomAD_AF_afr",
            "gnomAD_AF_amr",
            "gnomAD_AF_AsJ",
            "gnomAD_AF_nfe",
            "LMM_FLAGGED",
            "OMIM_GENESYMBOL",
            "MIMNUMBER",
            "OMIM_PHENO",
            "OR_V",
            "PHYLOP_PLCNTL",
            "PHYLOP_VRT",
            "GO_BP",
            "GO_CC",
            "GO_MF",
            "SP_BL",
            "SP_LODSCORE",
            "SVTYPE"
    );

    private final List<String> SNPSIFT_INFO = Arrays.asList(
            "VEP_canonical",
            "SIFT_scr",
            "SIFT_cnvtd_rnkscr",
            "SIFT_pred",
            "SIFT4G_scr",
            "SIFT4G_cnvtd_rnkscr",
            "SIFT4G_pred",
            "Polyphen2_HDIV_scr",
            "Polyphen2_HDIV_rnkscr",
            "Polyphen2_HDIV_pred",
            "Polyphen2_HVAR_scr",
            "Polyphen2_HVAR_rnkscr",
            "Polyphen2_HVAR_pred",
            "LRT_scr",
            "LRT_cnvtd_rnkscr",
            "LRT_pred",
            "LRT_Omega",
            "MutationTaster_scr",
            "MutationTaster_cnvtd_rnkscr",
            "MutationTaster_pred",
            "MutationTaster_model",
            "MutationTaster_AAE",
            "MutationAssessor_scr",
            "MutationAssessor_rnkscr",
            "MutationAssessor_pred",
            "FATHMM_scr",
            "FATHMM_cnvtd_rnkscr",
            "FATHMM_pred",
            "PROVEAN_scr",
            "PROVEAN_cnvtd_rnkscr",
            "PROVEAN_pred",
            "VEST4_scr",
            "VEST4_rnkscr",
            "MetaSVM_scr",
            "MetaSVM_rnkscr",
            "MetaSVM_pred",
            "MetaLR_scr",
            "MetaLR_rnkscr",
            "MetaLR_pred",
            "Reliability_index",
            "M-CAP_scr",
            "M-CAP_rnkscr",
            "M-CAP_pred",
            "REVEL_scr",
            "REVEL_rnkscr",
            "MutPred_scr",
            "MutPred_rnkscr",
            "MutPred_protID",
            "MutPred_AAchange",
            "MutPred_Top5features",
            "MVP_scr",
            "MVP_rnkscr",
            "MPC_scr",
            "MPC_rnkscr",
            "PrimateAI_scr",
            "PrimateAI_rnkscr",
            "PrimateAI_pred",
            "DEOGEN2_scr",
            "DEOGEN2_rnkscr",
            "DEOGEN2_pred",
            "BayesDel_addAF_scr",
            "BayesDel_addAF_rnkscr",
            "BayesDel_addAF_pred",
            "BayesDel_noAF_scr",
            "BayesDel_noAF_rnkscr",
            "BayesDel_noAF_pred",
            "LIST-S2_scr",
            "LIST-S2_rnkscr",
            "LIST-S2_pred",
            "Aloft_Fraction_transcripts_affected",
            "Aloft_prob_Tolerant",
            "Aloft_prob_Rec",
            "Aloft_prob_Dom",
            "Aloft_pred",
            "Aloft_Confidence",
            "DANN_scr",
            "DANN_rnkscr",
            "fathmm-MKL_coding_scr",
            "fathmm-MKL_coding_rnkscr",
            "fathmm-MKL_coding_pred",
            "fathmm-MKL_coding_group",
            "fathmm-XF_coding_scr",
            "fathmm-XF_coding_rnkscr",
            "fathmm-XF_coding_pred",
            "Eigen-raw_coding",
            "Eigen-raw_coding_rnkscr",
            "Eigen-phred_coding",
            "Eigen-PC-raw_coding",
            "Eigen-PC-raw_coding_rnkscr",
            "Eigen-PC-phred_coding",
            "GenoCanyon_scr",
            "GenoCanyon_rnkscr",
            "integrated_fitCons_scr",
            "integrated_fitCons_rnkscr",
            "integrated_confidence_value",
            "GM12878_fitCons_scr",
            "GM12878_fitCons_rnkscr",
            "GM12878_confidence_value",
            "H1-hESC_fitCons_scr",
            "H1-hESC_fitCons_rnkscr",
            "H1-hESC_confidence_value",
            "HUVEC_fitCons_scr",
            "HUVEC_fitCons_rnkscr",
            "HUVEC_confidence_value",
            "LINSIGHT",
            "LINSIGHT_rnkscr",
            "GERP___NR",
            "GERP___RS",
            "GERP___RS_rnkscr",
            "phyloP100way_vertebrate",
            "phyloP100way_vertebrate_rnkscr",
            "phyloP30way_mammalian",
            "phyloP30way_mammalian_rnkscr",
            "phyloP17way_primate",
            "phyloP17way_primate_rnkscr",
            "phastCons100way_vertebrate",
            "phastCons100way_vertebrate_rnkscr",
            "phastCons30way_mammalian",
            "phastCons30way_mammalian_rnkscr",
            "phastCons17way_primate",
            "phastCons17way_primate_rnkscr"
        );

    private final List<String> CASSENDRA_INFO = Arrays.asList(
            "RFG",
            "GI",
            "UCG",
            "ENS",
            "ENSGN",
            "MT",
            "MTGN",
            "CADD_RS",
            "CADD_PH",
            "PP_VB",
            "PP_PR",
            "PP_PL",
            "PC_VB",
            "PC_PR",
            "PC_PL",
            "MB",
            "FS_EN",
            "FS_SN",
            "FS_US",
            "FS_TG",
            "FS_SC",
            "FS_NS",
            "FS_WS",
            "TMAF",
            "TAMR",
            "TASN",
            "TAFR",
            "TEUR",
            "RSID",
            "CD",
            "NA",
            "NC",
            "NE",
            "NF",
            "NG",
            "NH",
            "NI",
            "NJ",
            "NK",
            "NL",
            "NM",
            "MS",
            "IET",
            "IEO",
            "IEN",
            "CG",
            "AV",
            "AT",
            "EV",
            "ET",
            "ENN",
            "ENC",
            "CID",
            "CGENE",
            "CSTRAND",
            "CCDS",
            "CAA",
            "CCNT",
            "ARIC_AA",
            "ARIC_EA",
            "ERBCTA_CT",
            "ERBCTA_NM",
            "ERBCTA_SC",
            "ERBSEG_CT",
            "ERBSEG_NM",
            "ERBSEG_SC",
            "ERBTFBS_TF",
            "ERBTFBS_PB",
            "SCSNV_ADA",
            "SCSNV_RS",
            "ERBSUM_NM",
            "ERBSUM_SC",
            "FE",
            "FC",
            "ENCSEG_CT",
            "ENCSEG_NM",
            "ENCDNA_SC",
            "ENCDNA_CT",
            "ENCTFBS_TF",
            "ENCTFBS_SC",
            "ENCTFBS_CL",
            "GRASP_RS",
            "GRASP_PMID",
            "GRASP_P",
            "GRASP_PH",
            "GRASP_AN",
            "GRASP_PL",
            "OREGANNO_TYPE",
            "OREGANNO_PMID",
            "RDB_MF",
            "RDB_WS",
            "SP_SC",
            "EX_AC",
            "EX_MQ",
            "GC",
            "GN",
            "SF",
            "SD",
            "SM",
            "SX",
            "OMIMS",
            "OMIMN",
            "OMIMT",
            "OMIMM",
            "OMIMC",
            "OMIMD",
            "OMIMMUS",
            "LiftedContig",
            "LiftedStart",
            "LiftedStop"
    );

    private long clinvar = 0L;
    private long cassandra = 0L;
    private long rejectedLiftover = 0L;
    private long funcotator = 0L;
    private long snpSift = 0L;

    private final VCFInfoHeaderLine UNABLE_TO_LIFT = new VCFInfoHeaderLine("LF", 1, VCFHeaderLineType.String, "Could not be lifted to alternate genome");

    private final Impact IMPACT_ANNOTATION = new Impact();

    private final Collection<String> ALLOWABLE_FILTERS = Arrays.asList("ReverseComplementedIndel", "NoTarget", "MismatchedRefAllele", "IndelStraddlesMultipleIntevals");

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);
        writer = createVCFWriter(new File(outFile));

        VCFHeader header = new VCFHeader(getHeaderForVariants());
        IMPACT_ANNOTATION.getDescriptions().forEach(header::addMetaDataLine);

        if (clinvarVariants != null) {
            VCFHeader clinvarHeader = (VCFHeader) getHeaderForFeatures(clinvarVariants);
            for (String id : CLINVAR_INFO) {
                VCFInfoHeaderLine line = clinvarHeader.getInfoHeaderLine(id);
                if (line == null) {
                    throw new GATKException("Clinvar missing expected header line: " + id);
                }
                header.addMetaDataLine(line);
            }

            List<String> allKeys = new ArrayList<>(clinvarHeader.getInfoHeaderLines().stream().map(VCFInfoHeaderLine::getID).toList());
            allKeys.removeAll(CLINVAR_INFO);
            if (!allKeys.isEmpty()) {
                logger.info("The following ClinVar fields will not be retained: " + StringUtils.join(allKeys, ", "));
            }
        }

        if (cassandraVariants != null) {
            VCFHeader cassandraHeader = (VCFHeader) getHeaderForFeatures(cassandraVariants);
            for (String id : CASSENDRA_INFO) {
                VCFInfoHeaderLine line = cassandraHeader.getInfoHeaderLine(id);
                if (line == null) {
                    logger.warn("Cassandra missing expected header line: " + id);
                    continue;
                }
                header.addMetaDataLine(line);
            }

            List<String> allKeys = new ArrayList<>(cassandraHeader.getInfoHeaderLines().stream().map(VCFInfoHeaderLine::getID).toList());
            allKeys.removeAll(CASSENDRA_INFO);
            if (!allKeys.isEmpty()) {
                logger.info("The following Cassandra fields will not be retained: " + StringUtils.join(allKeys, ", "));
            }
        }

        if (snpSiftVariants != null) {
            VCFHeader snpSiftHeader = (VCFHeader) getHeaderForFeatures(snpSiftVariants);
            for (String id : SNPSIFT_INFO) {
                VCFInfoHeaderLine line = snpSiftHeader.getInfoHeaderLine(id);
                if (line == null) {
                    logger.warn("SnpSift missing expected header line: " + id);
                    continue;
                }
                header.addMetaDataLine(line);
            }

            List<String> allKeys = new ArrayList<>(snpSiftHeader.getInfoHeaderLines().stream().map(VCFInfoHeaderLine::getID).toList());
            allKeys.removeAll(SNPSIFT_INFO);
            if (!allKeys.isEmpty()) {
                logger.info("The following SnpSift fields will not be retained: " + StringUtils.join(allKeys, ", "));
            }
        }

        if (funcotatorVariants != null) {
            VCFHeader funcotatorHeader = (VCFHeader) getHeaderForFeatures(funcotatorVariants);
            for (String id : FUNCOTATOR_INFO) {
                VCFInfoHeaderLine line = funcotatorHeader.getInfoHeaderLine(id);
                if (line == null) {
                    logger.warn("Funcotator missing expected header line: " + id);
                    continue;
                }
                header.addMetaDataLine(line);
            }

            List<String> allKeys = new ArrayList<>(funcotatorHeader.getInfoHeaderLines().stream().map(VCFInfoHeaderLine::getID).toList());
            allKeys.removeAll(FUNCOTATOR_INFO);
            if (!allKeys.isEmpty()) {
                logger.info("The following Funcotator fields will not be retained: " + StringUtils.join(allKeys, ", "));
            }
        }

        header.addMetaDataLine(UNABLE_TO_LIFT);

        writer.writeHeader(header);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);

        if (clinvarVariants != null) {
            for (VariantContext vc : featureContext.getValues(clinvarVariants)) {
                if (!matches(variant, vc)) {
                    continue;
                }

                clinvar++;
                transferInfoData(vcb, vc, CLINVAR_INFO);
            }
        }

        if (cassandraVariants != null) {
            for (VariantContext vc : featureContext.getValues(cassandraVariants)) {
                if (!matches(variant, vc)) {
                    continue;
                }

                cassandra++;
                transferInfoData(vcb, vc, CASSENDRA_INFO);
            }
        }

        if (snpSiftVariants != null) {
            for (VariantContext vc : featureContext.getValues(snpSiftVariants)) {
                if (!matches(variant, vc)) {
                    continue;
                }

                snpSift++;
                transferInfoData(vcb, vc, SNPSIFT_INFO);
            }
        }

        if (funcotatorVariants != null) {
            for (VariantContext vc : featureContext.getValues(funcotatorVariants)) {
                if (!matches(variant, vc)) {
                    continue;
                }

                funcotator++;
                transferInfoData(vcb, vc, FUNCOTATOR_INFO);
            }
        }

        if (liftoverRejectVariants != null) {
            for (VariantContext vc : featureContext.getValues(liftoverRejectVariants)) {
                if (!matches(variant, vc)) {
                    continue;
                }

                rejectedLiftover++;
                Set<String> filters = new HashSet<>(vc.getFilters());
                filters.retainAll(ALLOWABLE_FILTERS);
                vcb.attribute(UNABLE_TO_LIFT.getID(), StringUtils.join(filters, ","));
            }
        }

        Map<String, Object> toAnnotate = IMPACT_ANNOTATION.annotate(referenceContext, variant, null);
        if (toAnnotate != null) {
            vcb.putAttributes(toAnnotate);
        }

        writer.add(vcb.make());
    }

    private boolean matches(VariantContext source, VariantContext annotation){
        if (!source.getContig().equals(annotation.getContig()) &&
                source.getStart() == annotation.getStart() &&
                source.getEnd() == annotation.getEnd() &&
                source.getReference().equals(annotation.getReference())){
            logger.info("not matching: " + source.getContig() + "/" + source.getStart());
            return false;
        }

        //TODO: consider allele-specific annotations
        //if (!source.getAlleles().equals(annotation.getAlleles())){
        //    logger.warn("alleles do not match: " + source.getContig() + "/" + source.getStart());
        //}

        return true;
    }

    private void transferInfoData(VariantContextBuilder vcb, VariantContext source, List<String> ids){
        for (String id : ids){
            if (source.hasAttribute(id) && source.getAttribute(id) != null){
                vcb.attribute(id, source.getAttribute(id));
            }
        }
    }

    @Override
    public void closeTool(){
        writer.close();

        logger.info("Total variants annotated from ClinVar: " + clinvar);
        logger.info("Total variants annotated from Cassandra: " + cassandra);
        logger.info("Total variants annotated from Funcotator: " + funcotator);
        logger.info("Total variants annotated from SnpSift: " + snpSift);
        logger.info("Total variants annotated as failing liftover: " + rejectedLiftover);
    }
}
