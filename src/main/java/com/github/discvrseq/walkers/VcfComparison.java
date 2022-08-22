package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.ICSVWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This tool will compare the sites and genotypes between two VCFs. It will output site-level differences, including number of sites filtered/not between the query and reference VCF, and sites present/absent in either.
 * For each site that is present in both and not filtered, it will output the total number of genotypes that are discordant (which does not consider no-calls/filtered genotypes).
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar VcfComparison \
 *     -R currentGenome.fasta \
 *     -V myVCF.vcf \
 *     -rv referenceData.vcf.gz \
 *     -O outputTable.txt \
 *     -sites-output sitesTable.txt.gz
 *     --novel-or-altered-sites-vcf novelVariants.vcf.gz
 *     --missing-sites-vcf missingVariants.vcf.gz
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool will compare two VCFs at the site and genotype level, producing a summary table of concordance.",
        oneLineSummary = "Produces a concordance report between two VCFs at the site and genotype level",
        programGroup = VariantManipulationProgramGroup.class
)
public class VcfComparison extends ExtendedMultiVariantWalkerGroupedOnStart {
    @Argument(doc = "File to which the output table should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath outFile = null;

    @Argument(doc = "If provided, a tab-delimited list of each site with a discrepancy will be written", fullName = "sites-output", optional = true)
    public GATKPath siteOutput = null;

    @Argument(doc = "By default, any site from a VCF with genotypes, but where none of the samples are called, are treated as missing. Use this flag to include those sites", fullName = "include-no-call-sites", optional = true)
    public boolean includeNoCallSites = false;

    private ICSVWriter siteWriter = null;
    private VariantContextWriter novelSitesWriter;
    private VariantContextWriter missingSitesWriter;

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);
        IOUtil.assertFileIsWritable(outFile.toPath().toFile());

        if (siteOutput != null) {
            IOUtil.assertFileIsWritable(siteOutput.toPath().toFile());

            siteWriter = CsvUtils.getTsvWriter(siteOutput.toPath());
        }

        if (getVcfComparisonArgumentCollection().drivingVariantPaths.size() > 1) {
            throw new UserException.BadInput("This tool currently only supports one VCF input");
        }

        novelSitesWriter = initializeVcfWriter(getVcfComparisonArgumentCollection().novelSitesVcf, getDrivingVariantsFeatureInputs().get(0));
        missingSitesWriter = initializeVcfWriter(getVcfComparisonArgumentCollection().missingSitesVcf, getVcfComparisonArgumentCollection().refVariants);
    }

    private static final String INCLUSION_REASON = "IR";

    private VariantContextWriter initializeVcfWriter (GATKPath path, FeatureInput<VariantContext> vcfForHeader) {
        if (path == null) {
            return null;
        }

        IOUtil.assertFileIsWritable(path.toPath().toFile());

        VariantContextWriter writer = createVCFWriter(path);

        // setup the header fields
        VCFHeader header = (VCFHeader) getHeaderForFeatures(vcfForHeader);
        final Set<VCFHeaderLine> hInfo = new LinkedHashSet<>(header.getMetaDataInInputOrder());
        hInfo.add(new VCFInfoHeaderLine(INCLUSION_REASON, 1, VCFHeaderLineType.String, "The status of this variant relative to the reference VCF"));

        writer.writeHeader(new VCFHeader(hInfo, header.getGenotypeSamples()));

        return writer;
    }

    @Override
    protected VcfComparisonArgumentCollection getMultiVariantInputArgumentCollection() {
        return new VcfComparisonArgumentCollection();
    }

    private VcfComparisonArgumentCollection getVcfComparisonArgumentCollection() {
        return (VcfComparisonArgumentCollection)multiVariantInputArgumentCollection;
    }

    private static final class VcfComparisonArgumentCollection extends MultiVariantInputArgumentCollection.DefaultMultiVariantInputArgumentCollection {
        private static final long serialVersionUID = 1L;

        @Argument(doc = "Reference VCF", fullName = "referenceVCF", shortName = "rv", optional = false)
        public FeatureInput<VariantContext> refVariants = null;

        @Argument(doc = "Novel or Altered Sites VCF", fullName = "novel-or-altered-sites-vcf", shortName = "ns", optional = true)
        public GATKPath novelSitesVcf = null;

        @Argument(doc = "Missing Sites VCF", fullName = "missing-sites-vcf", shortName = "ms", optional = true)
        public GATKPath missingSitesVcf = null;

        @Override
        public List<GATKPath> getDrivingVariantPaths() {
            List<GATKPath> ret = new ArrayList<>(super.getDrivingVariantPaths());
            ret.add(refVariants);

            return ret;
        }
    }

    private int novelSitesRelativeToRef = 0;
    private int siteMissingRelativeToRef = 0;

    private int refFilteredNotSample = 0;
    private int sampleFilteredNotRef = 0;
    private int discordantGenotypes = 0;

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        Map<FeatureInput<VariantContext>, List<VariantContext>> variants = groupVariantsByFeatureInput(variantContexts);

        List<VariantContext> sampleVariants = variants.get(getDrivingVariantsFeatureInputs().get(0)).stream().filter(variantContext -> includeNoCallSites ||!variantContext.hasGenotypes() || variantContext.getCalledChrCount() > 0).collect(Collectors.toList());
        List<VariantContext> refVariants = variants.get(getVcfComparisonArgumentCollection().refVariants).stream().filter(variantContext -> includeNoCallSites || !variantContext.hasGenotypes() || variantContext.getCalledChrCount() > 0).collect(Collectors.toList());
        if (sampleVariants.isEmpty()) {
            possiblyWriteVariant(refVariants, referenceContext, "SiteMissingRelativeToRef", missingSitesWriter);
            siteMissingRelativeToRef++;
            return;
        }

        if (refVariants.isEmpty()){
            possiblyWriteVariant(sampleVariants, referenceContext, "NovelSitesRelativeToRef", novelSitesWriter);
            novelSitesRelativeToRef++;
            return;
        }


        if (sampleVariants.size() > 1) {
            logger.warn("More than one sample VCF record found for site: " + sampleVariants.get(0).getContig() + " " + sampleVariants.get(0).getStart());
        }

        if (refVariants.size() > 1) {
            logger.warn("More than one reference VCF record found for site: " + refVariants.get(0).getContig() + " " + refVariants.get(0).getStart());
        }

        boolean refIsFiltered = refVariants.stream().filter(variantContext -> !variantContext.isFiltered()).count() == 0;
        boolean sampleIsFiltered = refVariants.stream().filter(variantContext -> !variantContext.isFiltered()).count() == 0;

        if (sampleIsFiltered && refIsFiltered) {
            return;
        }
        else if (sampleIsFiltered && !refIsFiltered) {
            sampleFilteredNotRef++;
            possiblyWriteVariant(refVariants, referenceContext, "SampleFilteredNotRef", missingSitesWriter);
            return;
        }
        else if (!sampleIsFiltered && refIsFiltered) {
            refFilteredNotSample++;
            possiblyWriteVariant(sampleVariants, referenceContext, "RefFilteredNotSample", novelSitesWriter);
            return;
        }

        for (VariantContext vc : sampleVariants) {
            if (vc.isFiltered()) {
                continue;
            }

            for (String sn : vc.getSampleNames()) {
                Genotype g = vc.getGenotype(sn);
                if (g.isFiltered() || g.isNoCall()) {
                    continue;
                }

                for (VariantContext refVc : refVariants) {
                    if (refVc.isFiltered() || !refVc.hasGenotype(sn)) {
                        continue;
                    }

                    Genotype rg = refVc.getGenotype(sn);
                    if (rg.isFiltered() || rg.isNoCall()) {
                        continue;
                    }

                    if (!rg.sameGenotype(g)) {
                        possiblyWriteVariant(Collections.singletonList(vc), referenceContext, "Discordant Genotype: " + g.getSampleName(), novelSitesWriter);
                        discordantGenotypes++;
                    }
                }
            }
        }
    }

    private void possiblyWriteVariant(List<VariantContext> variants, ReferenceContext referenceContext, String message, @Nullable VariantContextWriter variantContextWriter) {
        if (siteWriter != null) {
            siteWriter.writeNext(new String[]{referenceContext.getContig(), String.valueOf(referenceContext.getStart()), message});
        }

        if (variantContextWriter != null) {
            variants.forEach(vc -> {
                VariantContextBuilder vcb = new VariantContextBuilder(vc);
                vcb.attribute(INCLUSION_REASON, message);

                variantContextWriter.add(vcb.make());
            });
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try (ICSVWriter writer = CsvUtils.getTsvWriter(outFile.toPath())){
            writer.writeNext(new String[]{"Metric", "Value"});

            writer.writeNext(new String[]{"NovelSitesRelativeToRef", String.valueOf(novelSitesRelativeToRef)});
            writer.writeNext(new String[]{"SitesMissingRelativeToRef", String.valueOf(siteMissingRelativeToRef)});
            writer.writeNext(new String[]{"RefFilteredNotSample", String.valueOf(refFilteredNotSample)});
            writer.writeNext(new String[]{"SampleFilteredNotRef", String.valueOf(sampleFilteredNotRef)});
            writer.writeNext(new String[]{"DiscordantGenotypes", String.valueOf(discordantGenotypes)});
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        try {
            try {
                if (siteWriter != null) {
                    siteWriter.close();
                }
            }
            catch (IOException e) {
                logger.error("Unable to close CSV writer", e);
            }

            try {
                if (novelSitesWriter != null) {
                    novelSitesWriter.close();
                }
            }
            catch (Exception e) {
                logger.error("Unable to close novel sites VCF writer", e);
            }

            try {
                if (missingSitesWriter != null) {
                    missingSitesWriter.close();
                }
            }
            catch (Exception e) {
                logger.error("Unable to close missing sites VCF writer", e);
            }
        }
        finally {
            super.closeTool();
        }
    }
}
