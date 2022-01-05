package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.ICSVWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

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
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool will compare two VCFs at the site and genotype level, producing a summary table of concordance.",
        oneLineSummary = "Produces a concordance report between two VCFs at the site and genotype level",
        programGroup = VariantManipulationProgramGroup.class
)
public class VcfComparison extends MultiVariantWalkerGroupedOnStart {
    @Argument(doc = "File to which the output table should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath outFile = null;

    @Argument(doc = "If provided, a tab-delimited list of each site with a discrepancy will be written", fullName = "sites-output", optional = true)
    public GATKPath siteOutput = null;

    private ICSVWriter siteWriter = null;

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

        // Cache map of source name -> FeatureInput
        drivingVariantSourceMap = new HashMap<>();
        getDrivingVariantsFeatureInputs().forEach(x -> drivingVariantSourceMap.put(x.getName(), x));
    }

    // maintain the mapping of source name (from VC) to FeatureInput name
    private Map<String, FeatureInput<VariantContext>> drivingVariantSourceMap;

    private Map<FeatureInput<VariantContext>, List<VariantContext>> groupVariantsByFeatureInput(final List<VariantContext> variants) {
        final Map<FeatureInput<VariantContext>, List<VariantContext>> byFeatureInput = new HashMap<>();
        variants.forEach(vc -> byFeatureInput.compute(drivingVariantSourceMap.get(vc.getSource()),
                (k, v) -> {
                    final List<VariantContext> variantList = v == null ? new ArrayList<>() : v;
                    variantList.add(vc);
                    return variantList;
                }
        ));

        for (FeatureInput<VariantContext> fi : getDrivingVariantsFeatureInputs()) {
            if (!byFeatureInput.containsKey(fi)) {
                byFeatureInput.put(fi, Collections.emptyList());
            }
        }
        return byFeatureInput;
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

        List<VariantContext> sampleVariants = variants.get(getDrivingVariantsFeatureInputs().get(0)).stream().filter(variantContext -> variantContext.getCalledChrCount() > 0).collect(Collectors.toList());
        if (sampleVariants.isEmpty()) {
            possiblyWriteVariant(referenceContext, "SiteMissingRelativeToRef");
            siteMissingRelativeToRef++;
            return;
        }

        List<VariantContext> refVariants = variants.get(getVcfComparisonArgumentCollection().refVariants).stream().filter(variantContext -> variantContext.getCalledChrCount() > 0).collect(Collectors.toList());
        if (refVariants.isEmpty()){
            possiblyWriteVariant(referenceContext, "NovelSitesRelativeToRef");
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
        boolean sampleIfFiltered = refVariants.stream().filter(variantContext -> !variantContext.isFiltered()).count() == 0;

        if (sampleIfFiltered && refIsFiltered) {
            return;
        }
        else if (sampleIfFiltered && !refIsFiltered) {
            sampleFilteredNotRef++;
            possiblyWriteVariant(referenceContext, "SampleFilteredNotRef");
            return;
        }
        else if (!sampleIfFiltered && refIsFiltered) {
            refFilteredNotSample++;
            possiblyWriteVariant(referenceContext, "RefFilteredNotSample");
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
                        possiblyWriteVariant(referenceContext, "Discordant Genotype: " + g.getSampleName());
                        discordantGenotypes++;
                    }
                }
            }
        }
    }

    private void possiblyWriteVariant(ReferenceContext referenceContext, String message) {
        if (siteWriter == null) {
            return;
        }

        siteWriter.writeNext(new String[]{referenceContext.getContig(), String.valueOf(referenceContext.getStart()), message});
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
            if (siteWriter != null) {
                siteWriter.close();
            }
        }
        catch (IOException e) {
            logger.error("Unable to close CSV writer", e);
        }
        finally {
            super.closeTool();
        }
    }
}
