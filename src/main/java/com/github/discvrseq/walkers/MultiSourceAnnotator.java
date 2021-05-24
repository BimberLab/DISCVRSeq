package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
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
        }

        //TODO: annotation source versions?
        //##META='Cassandra_version=15.4.10'
        //cassandraHeader.getMetaDataLine()

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
        logger.info("Total variants annotated as failing liftover: " + rejectedLiftover);
    }
}
