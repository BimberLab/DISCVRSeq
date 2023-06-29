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
            "dbNSFP_VEP_canonical",
            "dbNSFP_SIFT_scr",
            "dbNSFP_SIFT_cnvtd_rnkscr",
            "dbNSFP_SIFT_pred",
            "dbNSFP_SIFT4G_scr",
            "dbNSFP_SIFT4G_cnvtd_rnkscr",
            "dbNSFP_SIFT4G_pred",
            "dbNSFP_Polyphen2_HDIV_scr",
            "dbNSFP_Polyphen2_HDIV_rnkscr",
            "dbNSFP_Polyphen2_HDIV_pred",
            "dbNSFP_Polyphen2_HVAR_scr",
            "dbNSFP_Polyphen2_HVAR_rnkscr",
            "dbNSFP_Polyphen2_HVAR_pred",
            "dbNSFP_LRT_scr",
            "dbNSFP_LRT_cnvtd_rnkscr",
            "dbNSFP_LRT_pred",
            "dbNSFP_LRT_Omega",
            "dbNSFP_MutationTaster_scr",
            "dbNSFP_MutationTaster_cnvtd_rnkscr",
            "dbNSFP_MutationTaster_pred",
            "dbNSFP_MutationTaster_model",
            "dbNSFP_MutationTaster_AAE",
            "dbNSFP_MutationAssessor_scr",
            "dbNSFP_MutationAssessor_rnkscr",
            "dbNSFP_MutationAssessor_pred",
            "dbNSFP_FATHMM_scr",
            "dbNSFP_FATHMM_cnvtd_rnkscr",
            "dbNSFP_FATHMM_pred",
            "dbNSFP_PROVEAN_scr",
            "dbNSFP_PROVEAN_cnvtd_rnkscr",
            "dbNSFP_PROVEAN_pred",
            "dbNSFP_VEST4_scr",
            "dbNSFP_VEST4_rnkscr",
            "dbNSFP_MetaSVM_scr",
            "dbNSFP_MetaSVM_rnkscr",
            "dbNSFP_MetaSVM_pred",
            "dbNSFP_MetaLR_scr",
            "dbNSFP_MetaLR_rnkscr",
            "dbNSFP_MetaLR_pred",
            "dbNSFP_Reliability_index",
            "dbNSFP_M-CAP_scr",
            "dbNSFP_M-CAP_rnkscr",
            "dbNSFP_M-CAP_pred",
            "dbNSFP_REVEL_scr",
            "dbNSFP_REVEL_rnkscr",
            "dbNSFP_MutPred_scr",
            "dbNSFP_MutPred_rnkscr",
            "dbNSFP_MutPred_protID",
            "dbNSFP_MutPred_AAchange",
            "dbNSFP_MutPred_Top5features",
            "dbNSFP_MVP_scr",
            "dbNSFP_MVP_rnkscr",
            "dbNSFP_MPC_scr",
            "dbNSFP_MPC_rnkscr",
            "dbNSFP_PrimateAI_scr",
            "dbNSFP_PrimateAI_rnkscr",
            "dbNSFP_PrimateAI_pred",
            "dbNSFP_DEOGEN2_scr",
            "dbNSFP_DEOGEN2_rnkscr",
            "dbNSFP_DEOGEN2_pred",
            "dbNSFP_BayesDel_addAF_scr",
            "dbNSFP_BayesDel_addAF_rnkscr",
            "dbNSFP_BayesDel_addAF_pred",
            "dbNSFP_BayesDel_noAF_scr",
            "dbNSFP_BayesDel_noAF_rnkscr",
            "dbNSFP_BayesDel_noAF_pred",
            "dbNSFP_LIST-S2_scr",
            "dbNSFP_LIST-S2_rnkscr",
            "dbNSFP_LIST-S2_pred",
            "dbNSFP_Aloft_Fraction_transcripts_affected",
            "dbNSFP_Aloft_prob_Tolerant",
            "dbNSFP_Aloft_prob_Rec",
            "dbNSFP_Aloft_prob_Dom",
            "dbNSFP_Aloft_pred",
            "dbNSFP_Aloft_Confidence",
            "dbNSFP_DANN_scr",
            "dbNSFP_DANN_rnkscr",
            "dbNSFP_fathmm-MKL_coding_scr",
            "dbNSFP_fathmm-MKL_coding_rnkscr",
            "dbNSFP_fathmm-MKL_coding_pred",
            "dbNSFP_fathmm-MKL_coding_group",
            "dbNSFP_fathmm-XF_coding_scr",
            "dbNSFP_fathmm-XF_coding_rnkscr",
            "dbNSFP_fathmm-XF_coding_pred",
            "dbNSFP_Eigen-raw_coding",
            "dbNSFP_Eigen-raw_coding_rnkscr",
            "dbNSFP_Eigen-phred_coding",
            "dbNSFP_Eigen-PC-raw_coding",
            "dbNSFP_Eigen-PC-raw_coding_rnkscr",
            "dbNSFP_Eigen-PC-phred_coding",
            "dbNSFP_GenoCanyon_scr",
            "dbNSFP_GenoCanyon_rnkscr",
            "dbNSFP_integrated_fitCons_scr",
            "dbNSFP_integrated_fitCons_rnkscr",
            "dbNSFP_integrated_confidence_value",
            "dbNSFP_GM12878_fitCons_scr",
            "dbNSFP_GM12878_fitCons_rnkscr",
            "dbNSFP_GM12878_confidence_value",
            "dbNSFP_H1-hESC_fitCons_scr",
            "dbNSFP_H1-hESC_fitCons_rnkscr",
            "dbNSFP_H1-hESC_confidence_value",
            "dbNSFP_HUVEC_fitCons_scr",
            "dbNSFP_HUVEC_fitCons_rnkscr",
            "dbNSFP_HUVEC_confidence_value",
            "dbNSFP_LINSIGHT",
            "dbNSFP_LINSIGHT_rnkscr",
            "dbNSFP_GERP___NR",
            "dbNSFP_GERP___RS",
            "dbNSFP_GERP___RS_rnkscr",
            "dbNSFP_phyloP100way_vertebrate",
            "dbNSFP_phyloP100way_vertebrate_rnkscr",
            "dbNSFP_phyloP30way_mammalian",
            "dbNSFP_phyloP30way_mammalian_rnkscr",
            "dbNSFP_phyloP17way_primate",
            "dbNSFP_phyloP17way_primate_rnkscr",
            "dbNSFP_phastCons100way_vertebrate",
            "dbNSFP_phastCons100way_vertebrate_rnkscr",
            "dbNSFP_phastCons30way_mammalian",
            "dbNSFP_phastCons30way_mammalian_rnkscr",
            "dbNSFP_phastCons17way_primate",
            "dbNSFP_phastCons17way_primate_rnkscr"
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
