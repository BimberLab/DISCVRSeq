package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.ICSVWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
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

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;


/**
 * This tool will compare the samples in a VCF against and number of reference VCFs.  It will produce a report summarizing the number of alleles shared between each sample and the reference sets, based on AF
 * One possible use of this tool is to compare VCF samples against reference panels of alleles, such as alleles unique to different ethnicities or populations.
 *
 * Each of the reference VCF(s) must have an AF annotation, representing the frequency of each allele in that population. Per sample, the tool tracks the cumulative sum of AFs for each matching site/allele.
 * Per sample, it also tracks the min/max possible score that could be obtained based on overlapping sites, and the number of overlapping sites. This is used to produce a final 0-1 score.
 *
 *
 * <h3>Usage example (note naming of the --ref-sites files):</h3>
 * <pre>
 *  java -jar DISCVRseq.jar VariantConcordanceScore \
 *     -R reference.fasta \
 *     --ref-sites:SET1 ref1.vcf \
 *     --ref-sites:SET2 ref2.vcf.gz \
 *     -V myVCF.vcf.gz \
 *     -O output.txt
 * </pre>
 *
 * This will output a table with one row for each sample/reference pair, with the following data:
 * <ul>
 *  <li>SampleName: name of sample</li>
 *  <li>ReferenceName: name of the reference, as given on the command line ("i.e. --ref-sites:SET1 ref1.vcf)</li>
 *  <li>MarkersWithData: The number of sites from this sample with genotype data overlapping this reference</li>
 *  <li>MarkersNoData: The number of sites from this reference where the sample lacked genotype data</li>
 *  <li>FractionWithData: The fraction of sites from the reference where the sample had genotype data</li>
 *  <li>UniqueContigs: The number of unique contigs where a marker overlapped in this sample</li>
 *  <li>CumulativeAF: The cumulative allele frequency for the alles of all genotypes overlapping this reference</li>
 *  <li>MinPossible: The minimum possible AF score that could be achieved for these sites</li>
 *  <li>MaxPossible: The maximum possible AF score that could be achieved for these sites</li>
 *  <li>Score: The scaled score, which is: (CumulativeAF-MinPossible) / (MaxPossible-MinPossible)</li>
 * </ul>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool compares an input VCF against a set of reference VCFs, reporting concordance between them",
        oneLineSummary = "Summarize allele-level concordance between a VCF and reference VCFs",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class VariantConcordanceScore extends ExtendedMultiVariantWalkerGroupedOnStart {

    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsWritable(new File(outFile));

        for (FeatureInput<VariantContext> ref : getVariantConcordanceScoreArgumentCollection().referenceFiles)
        {
            VCFHeader header = (VCFHeader)getHeaderForFeatures(ref);
            if (!header.hasInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                throw new GATKException("VCF missing AF annotation: " + ref.getName());
            }

            if (!ref.hasUserSuppliedName()) {
                throw new GATKException("A unique name should be supplied for all references on the command line (i.e. '--ref-sites:SET1 ref1.vcf'), missing for: " + ref.getName());
            }
        }

        VCFHeader header = (VCFHeader)getHeaderForFeatures(getVariantConcordanceScoreArgumentCollection().inputVariants);
        header.getSampleNamesInOrder().forEach(s -> sampleMap.put(s, new ScoreBySample(s)));
    }

    private final Map<String, ScoreBySample> sampleMap = new HashMap<>();

    private static class ScoreBySample {
        final String id;
        final Map<String, ScoreByPopulation> populations = new HashMap<>();

        public ScoreBySample(String id) {
            this.id = id;
        }
    }

    private static class ScoreByPopulation {
        Double cumulativeAfScore = 0.0;
        Double minPossibleAfScore = 0.0;
        Double maxPossibleAfScore = 0.0;
        Long markersWithData = 0L;
        Long markersNoData = 0L;
        Set<String> uniqueContigs = new HashSet<>();
    }

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        Map<FeatureInput<VariantContext>, List<VariantContext>> grouped = groupVariantsByFeatureInput(variantContexts);

        boolean hasSites = false;
        for (FeatureInput<VariantContext> feat : getVariantConcordanceScoreArgumentCollection().referenceFiles) {
            if (grouped.get(feat).size() > 1) {
                throw new GATKException("Duplicate variants detected in reference VCF: " + feat.getName() + " for position: " + variantContexts.get(0).getContig() + ": " + variantContexts.get(0).getStart());
            }
            else if (!grouped.get(feat).isEmpty()) {
                hasSites = true;
            }
        }

        if (!hasSites) {
            return;
        }

        List<VariantContext> sampleVariants = grouped.get(getVariantConcordanceScoreArgumentCollection().inputVariants);
        sampleVariants = sampleVariants.stream().filter(vc -> !vc.isFiltered()).collect(Collectors.toList());
        if (sampleVariants.isEmpty()) {
            return;
        }
        else if (sampleVariants.size() > 1) {
            // There are rare instances where the VCF has two lines for the same site, with different REF alleles (i.e. C, and CTG).
            // This is probably from suboptimal merges between two VCFs. For now, just limit to the site that matches this REF
            logger.warn("Duplicate variants detected in input VCF: " + " for position: " + variantContexts.get(0).getContig() + ": " + variantContexts.get(0).getStart());
            return;
        }

        for (FeatureInput<VariantContext> population : getVariantConcordanceScoreArgumentCollection().referenceFiles) {
            if (!grouped.containsKey(population) || grouped.get(population).isEmpty()) {
                continue;
            }

            processPopulation(population, grouped.get(population).get(0), sampleVariants.get(0));
        }
    }

    private void processPopulation(FeatureInput<VariantContext> population, VariantContext referenceSite, VariantContext sampleSite) {
        Map<Allele, Double> afMap = new HashMap<>();
        final List<Double> afVals = referenceSite.getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0);
        if (afVals.size() != referenceSite.getAlternateAlleles().size())
        {
            throw new GATKException("AFs did not match alleles in ref: " + population.getName() + " at position: " + referenceSite.getStart());
        }

        referenceSite.getAlternateAlleles().forEach(a -> afMap.put(a, afVals.get(referenceSite.getAlleleIndex(a) - 1)));

        double refAF = 1.0;
        for (Double d : afVals) {
            refAF = refAF - d;
        }

        if (refAF < 0) {
            throw new GATKException("Improper AF for reference: " + population.getName() + ", at site: " + referenceSite.getContig() + ": " + referenceSite.getStart());
        }
        afMap.put(referenceSite.getReference(), refAF);
        afVals.add(refAF);

        final double maxScore = Collections.max(afVals);
        final double minScore = Collections.min(afVals);

        sampleSite.getSampleNames().forEach(sn -> {
            ScoreBySample s = sampleMap.get(sn);
            ScoreByPopulation sp = s.populations.containsKey(population.getName()) ? s.populations.get(population.getName()) : new ScoreByPopulation();
            Genotype g = sampleSite.getGenotype(sn);
            if (!g.isCalled() || g.isFiltered()) {
                sp.markersNoData++;

            }
            else {
                double minScoreForAnimal = 0.0;
                for (Allele a : g.getAlleles()) {
                    if (afMap.containsKey(a)) {
                        minScoreForAnimal += minScore;
                    }

                    sp.cumulativeAfScore += afMap.getOrDefault(a, 0.0);
                }

                sp.markersWithData++;
                sp.uniqueContigs.add(referenceSite.getContig());
                sp.minPossibleAfScore += minScoreForAnimal;
                sp.maxPossibleAfScore += maxScore * 2;
            }

            s.populations.put(population.getName(), sp);
        });
    }

    @Override
    public Object onTraversalSuccess() {
        //Note: need should consider implementing some kind of logic to make actual calls per reference
        NumberFormat format = NumberFormat.getInstance();
        format.setMaximumFractionDigits(3);

        try (ICSVWriter output = CsvUtils.getTsvWriter(new File(outFile))) {
            output.writeNext(new String[]{"SampleName", "ReferenceName", "MarkersWithData", "MarkersNoData", "FractionWithData", "UniqueContigs", "CumulativeAF", "MinPossible", "MaxPossible", "Score"});

            for (String sample : sampleMap.keySet()) {
                ScoreBySample ss = sampleMap.get(sample);
                for (FeatureInput<VariantContext> ref : getVariantConcordanceScoreArgumentCollection().referenceFiles) {
                    ScoreByPopulation sp = ss.populations.get(ref.getName());
                    if (sp == null) {
                        sp = new ScoreByPopulation();
                    }

                    double fractionWithData = sp.markersWithData == 0 ? 0.0 : (double)sp.markersWithData / (double)(sp.markersWithData + sp.markersNoData);
                    double score = sp.cumulativeAfScore == 0.0 ? 0.0 : (sp.cumulativeAfScore - sp.minPossibleAfScore) / (sp.maxPossibleAfScore - sp.minPossibleAfScore);
                    output.writeNext(new String[]{
                            sample, ref.getName(),
                            String.valueOf(sp.markersWithData),
                            String.valueOf(sp.markersNoData),
                            format.format(fractionWithData),
                            String.valueOf(sp.uniqueContigs.size()),
                            String.valueOf(sp.cumulativeAfScore),
                            String.valueOf(sp.minPossibleAfScore),
                            String.valueOf(sp.maxPossibleAfScore),
                            format.format(score)
                    });
                }
            }

        }
        catch (IOException e) {
            throw new GATKException(e.getMessage());
        }

        return super.getTraversalIntervals();
    }

    @Override
    protected VariantConcordanceScoreArgumentCollection getMultiVariantInputArgumentCollection() {
        return new VariantConcordanceScoreArgumentCollection();
    }

    private VariantConcordanceScoreArgumentCollection getVariantConcordanceScoreArgumentCollection() {
        return (VariantConcordanceScoreArgumentCollection)multiVariantInputArgumentCollection;
    }

    private static final class VariantConcordanceScoreArgumentCollection extends MultiVariantInputArgumentCollection {
        private static final long serialVersionUID = 1L;

        @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
                doc = "One or more VCF files containing variants", common = false, optional = false)
        public FeatureInput<VariantContext> inputVariants;

        @Argument(fullName = "ref-sites", shortName = "rs", doc = "VCF file containing sites to test.  Must be uniquely named", optional = false)
        public List<FeatureInput<VariantContext>> referenceFiles = new ArrayList<>();

        @Override
        public List<GATKPath> getDrivingVariantPaths() {
            List<GATKPath> ret = new ArrayList<>();
            ret.add(inputVariants);
            ret.addAll(referenceFiles);

            return ret;
        }
    }
}
