package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import com.github.discvrseq.walkers.annotator.Impact;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
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
            "FATHMM_NonCodingScore",
            "FATHMM_NonCodingGroup",
            "FATHMM_SCORE",
            "FATHMM_CodingGroup",
            "CADD_PHRED_Quality",
            "CAPICE_Score",
            "ClinPredScore",
            "TargetSiteDuplication",
            "SVTYPE",
            "MitoC_START",
            "Mito_LEN",
            "Mito_END",
            "MobileElementINFO",
            "RefAltDiff",
            "StracturalVariantionIMPRECISE",
            "CI_END",
            "CI_POS",
            "END",
            "CallSet",
            "MergedCalls",
            "SVOverlap_AF",
            "NumberSamples",
            "SVOverlap_EAS_AF",
            "SVOverlap_EUR_AF",
            "SVOverlap_AFR_AF",
            "SVOverlap_AMR_AF",
            "SVOverlap_SAS_AF",
            "SiteMedPosteriorProb",
            "ARIC_score",
            "DBSCSNV_splicingConsensus",
            "FANTOM_ENHANCERID",
            "FANTOM_SCORE",
            "STRAND",
            "DNASE1_hypersinsitivitysite",
            "GRASP_phenotypePvalue",
            "phylopPlacentalSCORE",
            "PHYLOPVertScore",
            "LODSCORE",
            "BranchLength",
            "RS_dbSNP_ID",
            "CHR_POS",
            "ReverseStrand",
            "VariationProperty",
            "GENEINFO",
            "dbSNPBUildID",
            "VariantAlleleOrigin",
            "SSR",
            "Weight",
            "VariantClass",
            "PubMed",
            "ThridPartyAnn",
            "PubMedLink",
            "Has3DStructure",
            "SubmitterLinkOut",
            "NonSynFramshift",
            "NonSynMissense",
            "NonSynNonsense",
            "REFhasIdentical",
            "SYN_codingVariation",
            "3UTR",
            "5UTR",
            "AcceptorSpliceSite",
            "DonorSpliceSite",
            "INTON",
            "3PrimeGeneRegion",
            "5PrimeGeneRegion",
            "AnnotherVariantMaps",
            "HasAssemblyConflict",
            "AssemblySpecific",
            "MUTATION",
            "Validated",
            "Greaterthan5AlleleFreqinAllPops",
            "Greaterthan5AlleleFreqinAPop",
            "MarkerOnGenotypingKit",
            "GenotypesAvailable",
            "1000GenomesPhase1",
            "1000GenomesPhase3",
            "ClinicalDiagnosticAssay",
            "LocusSpecificDB",
            "MicroAttThirdParty",
            "HasOmim",
            "ContigAlleleNotPresent",
            "WithdrawnBySubmitter",
            "NonOverlappingAllele",
            "NonCoding",
            "AFby1000Genomes",
            "CommonSNP",
            "TOPMED",
            "ENSEMBLEGENE",
            "ENSEMBLEID",
            "MIMNUMBER",
            "GENESYMBOL",
            "PHENOTYPES",
            "Top_Genes",
            "LOF_Mechanism",
            "Mode_of_Inheritance",
            "GenesNotesLOF",
            "ACMG59rec",
            "Symbol",
            "GeneID",
            "Chr",
            "Chr_Band",
            "Cancer_Somatic_Mut",
            "Cancer_Germline_Mut",
            "Tumour_Types_Somatic",
            "Tumour_Types_germline",
            "Cancer_Syndrome",
            "Tissue_Types",
            "Cancer_Molecular_Genetics",
            "Mutation _Type",
            "Translocation_Partner",
            "Other_Germline_Mut",
            "Other_Syndrome",
            "AF_from GOESP",
            "AF_EXAC",
            "AF_TGP",
            "ALLELEID",
            "CLNDiseaseName",
            "CLNDVARName",
            "CLNDVarPair",
            "CLNDVarPairDisease",
            "CLN_HGVS_expression",
            "CLNVarReviewStatus",
            "ClinicalSignifigance",
            "ClinicalSignifiganceConflict",
            "ClinicalSigHaplotype",
            "CLNV_VariantType",
            "CLNV_OntologyID",
            "CLNV_ID_pairs",
            "DBVAR_ID",
            "GENEINFO",
            "MolecularConsequence",
            "ORIGIN",
            "RS_DBSNP_ID",
            "SSR",
            "overlapping_mutations",
            "total_alterations_in_gene",
            "fusion_genes",
            "tissue_types_affected",
            "AF_gnomAD",
            "gnomAD_AF_afr",
            "gnomAD_AF_afr_female",
            "gnomAD_AF_afr_male",
            "gnomAD_AF_amr",
            "gnomAD_AF_amr_female",
            "gnomAD_AF_amr_male",
            "gnomAD_AF_AshkenaziJ",
            "gnomAD_AF_asj_female",
            "gnomAD_AS_asj_male",
            "gnomAD_AF_eas",
            "gnomAD_AF_eas_female",
            "gnomAD_AF_eas_male",
            "gnomAD_AF_female",
            "gnomAD_AF_fin",
            "gnomAD_AF_fin_female",
            "gnomAD_AF_fin_male",
            "gnomAD_AF_male",
            "gnomAD_AF_nfe",
            "gnomAD_AF_nfe_est",
            "gnomAD_AF_nfe_female",
            "gnomAD_AF_nfe_male",
            "gnomAD_AF_nfe_nwe",
            "gnomAD_AF_nfe_onf",
            "gnomAD_AF_nfe_seu",
            "gnomAD_AF_oth",
            "gnomAD_AF_oth_female",
            "gnomAD_AF_oth_male",
            "gnomAD_AF_pop_max",
            "gnomAD_AF_raw",
            "LMM_FLAGGED",
            "Reference",
            "Synonym Gene",
            "familial_Syndrome",
            "mRNA_id prot_acc",
            "ValuesOfKnownReg",
            "Uniprot_gene",
            "uniprot_entry_name",
            "DrugBank",
            "alt_uniprot_accessions",
            "uniprot_accession",
            "GO_Biolocial_Process",
            "GO_Cellular_component",
            "GO_Molecular_Function",
            "TranscritptionFactor",
            "Strand",
            "FANTOM_TFBS_SCORE",
            "Funseq2_Score",
            "EnsembleType"
    );

    private final List<String> SNPSIFT_INFO = Arrays.asList(
            "dbNSFP_1000Gp3_AC",
            "dbNSFP_1000Gp3_AF",
            "dbNSFP_1000Gp3_AFR_AC",
            "dbNSFP_1000Gp3_AFR_AF",
            "dbNSFP_1000Gp3_AMR_AC",
            "dbNSFP_1000Gp3_AMR_AF",
            "dbNSFP_1000Gp3_EAS_AC",
            "dbNSFP_1000Gp3_EAS_AF",
            "dbNSFP_1000Gp3_EUR_AC",
            "dbNSFP_1000Gp3_EUR_AF",
            "dbNSFP_1000Gp3_SAS_AC",
            "dbNSFP_1000Gp3_SAS_AF",
            "dbNSFP_CADD_phred",
            "dbNSFP_ESP6500_AA_AC",
            "dbNSFP_ESP6500_AA_AF",
            "dbNSFP_ESP6500_EA_AC",
            "dbNSFP_ESP6500_EA_AF",
            "dbNSFP_ExAC_AC",
            "dbNSFP_ExAC_AF",
            "dbNSFP_ExAC_AFR_AC",
            "dbNSFP_ExAC_AFR_AF",
            "dbNSFP_ExAC_AMR_AC",
            "dbNSFP_ExAC_AMR_AF",
            "dbNSFP_ExAC_Adj_AC",
            "dbNSFP_ExAC_Adj_AF",
            "dbNSFP_ExAC_EAS_AC",
            "dbNSFP_ExAC_EAS_AF",
            "dbNSFP_ExAC_FIN_AC",
            "dbNSFP_ExAC_FIN_AF",
            "dbNSFP_ExAC_NFE_AC",
            "dbNSFP_ExAC_NFE_AF",
            "dbNSFP_ExAC_SAS_AC",
            "dbNSFP_ExAC_SAS_AF",
            "dbNSFP_FATHMM_pred",
            "dbNSFP_GERP___NR",
            "dbNSFP_GERP___RS",
            "dbNSFP_Interpro_domain",
            "dbNSFP_LRT_pred",
            "dbNSFP_MetaSVM_pred",
            "dbNSFP_MutationAssessor_pred",
            "dbNSFP_MutationTaster_pred",
            "dbNSFP_PROVEAN_pred",
            "dbNSFP_Polyphen2_HDIV_pred",
            "dbNSFP_Polyphen2_HVAR_pred",
            "dbNSFP_SIFT_pred",
            "dbNSFP_Uniprot_acc",
            "dbNSFP_phastCons100way_vertebrate"
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

    private final Impact IMPACT = new Impact();

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

        Map<String, Object> toAnnotate = IMPACT.annotate(referenceContext, variant, null);
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
