package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFMetaHeaderLine;
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

    );

    private final List<String> SNPSIFT_INFO = Arrays.asList(
            "1000Gp1_AF",
            "1000Gp1_AFR_AF",
            "1000Gp1_AMR_AF",
            "1000Gp1_ASN_AF",
            "1000Gp1_EUR_AF",
            "1000Gp3_AC",
            "1000Gp3_AF",
            "1000Gp3_AFR_AC",
            "1000Gp3_AFR_AF",
            "1000Gp3_AMR_AC",
            "1000Gp3_AMR_AF",
            "1000Gp3_EAS_AC",
            "1000Gp3_EAS_AF",
            "1000Gp3_EUR_AC",
            "1000Gp3_EUR_AF",
            "1000Gp3_SAS_AC",
            "1000Gp3_SAS_AF",
            "CADD_phred",
            "ESP6500_AA_AC",
            "ESP6500_AA_AF",
            "ESP6500_AA_AF",
            "ESP6500_EA_AC",
            "ESP6500_EA_AF",
            "ESP6500_EA_AF",
            "ExAC_AC",
            "ExAC_Adj_AC",
            "ExAC_Adj_AF",
            "ExAC_AF",
            "ExAC_AFR_AC",
            "ExAC_AFR_AF",
            "ExAC_AMR_AC",
            "ExAC_AMR_AF",
            "ExAC_EAS_AC",
            "ExAC_EAS_AF",
            "ExAC_FIN_AC",
            "ExAC_FIN_AF",
            "ExAC_NFE_AC",
            "ExAC_NFE_AF",
            "ExAC_SAS_AC",
            "ExAC_SAS_AF",
            "FATHMM_pred",
            "GERP++_NR",
            "GERP++_RS",
            "Interpro_domain",
            "LRT_pred",
            "MetaSVM_pred",
            "MutationAssessor_pred",
            "MutationTaster_pred",
            "MutationTaster_pred",
            "phastCons100way_vertebrate",
            "Polyphen2_HDIV_pred",
            "Polyphen2_HVAR_pred",
            "PROVEAN_pred",
            "SIFT_pred",
            "Uniprot_acc"
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

    private final Collection<String> ALLOWABLE_FILTERS = Arrays.asList("ReverseComplementedIndel", "NoTarget", "MismatchedRefAllele", "IndelStraddlesMultipleIntevals");

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);
        writer = createVCFWriter(new File(outFile));

        VCFHeader header = new VCFHeader(getHeaderForVariants());

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
                    throw new GATKException("Cassandra missing expected header line: " + id);
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
                    throw new GATKException("SnpSift missing expected header line: " + id);
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
                    throw new GATKException("Funcotator missing expected header line: " + id);
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
