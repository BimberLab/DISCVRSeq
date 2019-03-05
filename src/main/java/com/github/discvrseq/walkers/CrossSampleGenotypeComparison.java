package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.opencsv.CSVWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This tool was originally created to help identify possible sample mix-ups. It accepts an input VCF, reference VCF and list of samples to consider.  It will compare each sample in this list from the VCF against every samples in the reference VCF, 
 * and produce a summary table about concordance.
 *
 * <h3>Usage example Without Sample List (will compare all):</h3>
 * <pre>
 *  java -jar DISCVRseq.jar CrossSampleGenotypeComparison \
 *     -R currentGenome.fasta \
 *     -V myVCF.vcf \
 *     -rv referenceData.vcf.gz \
 *     -O outputTable.txt \
 *     --summaryTable summaryTable.txt
 * </pre>
 * <h3>Usage Example With Sample List:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar CrossSampleGenotypeComparison \
 *     -R currentGenome.fasta \
 *     -V myVCF.vcf \
 *     -rv referenceData.vcf.gz \
 *     -O outputTable.txt \
 *     --summaryTable summaryTable.txt \
 *     -s Sample1 \
 *     -s Sample2 \
 *     -s Sample13
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool was originally created to help detect potential sample mix-up.  It accepts an input VCF, reference VCF and list of samples to consider.  It will compare each sample in this list from the VCF against every samples in the reference VCF, and produce a summary table about concordance.",
        oneLineSummary = "Produces a concordance report between samples in a VCF against every other sample in a reference VCF",
        programGroup = DiscvrSeqProgramGroup.class
)
public class CrossSampleGenotypeComparison extends VariantWalker {
    @Argument(doc = "File to which the output table should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(doc = "File to which the summary table should be written", fullName = "summaryTable", shortName = "st", optional = false)
    public String summaryOutFile = null;

    @Argument(doc = "Reference VCF", fullName = "referenceVCF", shortName = "rv", optional = false)
    public FeatureInput<VariantContext> refVariants = null;

    @Argument(doc = "Samples To Consider", fullName = "samples", shortName = "s", optional = true)
    public List<String> SAMPLES = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        if (SAMPLES.isEmpty()) {
            SAMPLES.addAll(getHeaderForVariants().getSampleNamesInOrder());
        }

        VCFHeader header = (VCFHeader)getHeaderForFeatures(refVariants);
        refSamples.addAll(header.getSampleNamesInOrder());
    }

    private List<String> refSamples = new ArrayList<>();
    private int warningsLogged = 0;
    private Map<String, Tracker> results = new HashMap<>();

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (vc.isFiltered()) {
            return;
        }

        //TODO: ignored list?
        if (vc.getContig().equalsIgnoreCase("chrUn")){
            return;
        }

        for (String sn : SAMPLES) {
            List<VariantContext> list = featureContext.getValues(refVariants);
            list.removeIf(x -> x.getStart() != vc.getStart());

            if (list.isEmpty()){
                if (warningsLogged < 10) {
                    logger.warn("position not found in reference VCF: " + vc.getContig() + ":" + vc.getStart());
                    warningsLogged++;

                    if (warningsLogged == 10){
                        logger.warn("future warnings will not be logged");
                    }
                }

                return;
            }

            if (vc.hasGenotype(sn)) {
                Genotype g = vc.getGenotype(sn);
                if (g.isFiltered() || g.isNoCall()) {
                    continue;
                }

                //iterate all animals in the ref:
                for (VariantContext c : list){
                    for (String sn2 : refSamples) {
                        String key = Tracker.getKey(sn, sn2);
                        Tracker t = results.get(key);
                        if (t == null){
                            t = new Tracker(sn, sn2);
                            results.put(key, t);
                        }

                        t.totalGenotypes++;
                        Genotype refGenotype = c.getGenotype(sn2);
                        if (refGenotype != null && !refGenotype.isFiltered() && !refGenotype.isNoCall()) {
                            if (!refGenotype.sameGenotype(g)) {
                                t.totalDiscordant++;
                            }
                            else {
                                t.totalConcordant++;
                            }
                        }
                    }
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(new File(outFile)), '\t', CSVWriter.NO_ESCAPE_CHARACTER, CSVWriter.NO_QUOTE_CHARACTER); CSVWriter summaryWriter = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(new File(summaryOutFile)), '\t', CSVWriter.NO_ESCAPE_CHARACTER, CSVWriter.NO_QUOTE_CHARACTER)){
            writer.writeNext(new String[]{"SampleName", "RefSampleName", "TotalCalledGenotypes", "TotalComparedToRef", "TotalDiscordant", "FractionDiscordant"});

            List<String> summaryHeader = new ArrayList<>();
            summaryHeader.add("");
            summaryHeader.addAll(refSamples);
            summaryWriter.writeNext(summaryHeader.toArray(new String[summaryHeader.size()]));

            for (String sn : SAMPLES){
                List<String> summaryLine = new ArrayList<>();
                summaryLine.add(sn);

                for (String sn2 : refSamples){
                    Tracker t = results.get(Tracker.getKey(sn, sn2));
                    if (t == null){
                        summaryLine.add("");
                        continue;
                    }

                    int totalCompared = t.totalConcordant + t.totalDiscordant;
                    String fraction = (totalCompared == 0 ? "" : String.valueOf(t.totalDiscordant / (double)totalCompared));
                    writer.writeNext(new String[]{t.sample1, t.sample2, String.valueOf(t.totalGenotypes), String.valueOf(t.totalConcordant + t.totalDiscordant), String.valueOf(t.totalDiscordant), fraction});

                    summaryLine.add(fraction);
                }

                summaryWriter.writeNext(summaryLine.toArray(new String[summaryLine.size()]));
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return super.onTraversalSuccess();
    }

    private static class Tracker {
        String sample1;
        String sample2;

        int totalGenotypes = 0;
        int totalConcordant = 0;
        int totalDiscordant = 0;

        private Tracker(String sample1, String sample2) {
            this.sample1 = sample1;
            this.sample2 = sample2;
        }

        private static String getKey(String sample1, String sample2) {
            return sample1 + "||" + sample2;
        }
    }
}
