package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;


/**
 * This walker will iterate a BAM and report sites with passing alignments with reads from more than 2 alleles, where
 * the third allele is above the supplied threshold.  The original purpose was to find positions of potential duplications
 * or poor genome quality (i.e. reads from more than one duplicated gene are mapping to a single site).
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar MultipleAllelesAtLoci \
 *     -R reference.fasta \
 *     -I bams.list \
 *     -O output.bed
 * </pre>
 * <p></p>
 * would produce a BED file that looks like:
 * <p></p>
 * <pre>
 *     chr20 10000000 10000001 MultiAllelicSite 2   +   Subj1: [A/0.25/1/4];Subj12: [A/0.2/2/10];
 *     chr20 10000865 10000866 MultiAllelicSite 1   +   Subj4: [T/0.1/1/10]
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool iterates a set of BAMs, reporting any loci (per sample) with more then three bases, based on the supplied threshold.",
        oneLineSummary = "Generate a list of loci with three or more distinct bases detected",
        programGroup = DiscvrSeqProgramGroup.class
)
public class MultipleAllelesAtLoci extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The output BED file", common = false, optional = true)
    private File outputFile = null;

    private PrintStream outputStream = null;

    /**
     * Reads with MAPQ > minMappingQuality are treated as usable for variation detection, contributing to the PASS
     * state.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth.", optional = true)
    byte minMappingQuality = 10;

    /**
     * Bases with less than minBaseQuality are viewed as not sufficiently high quality to contribute to the PASS state
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth.", optional = true)
    byte minBaseQuality = 20;

    /**
     * If the number of QC+ bases (on reads with MAPQ > minMappingQuality and with base quality > minBaseQuality) exceeds this
     * value and is less than maxDepth the site is considered PASS.
     */
    @Advanced
    @Argument(fullName = "minDepth", shortName = "minDepth", doc = "Minimum QC+ read depth before a locus is considered callable", optional = true)
    int minDepth = 4;

    @Argument(fullName = "minBasePct", shortName = "minBasePct", doc = "If a given site has more than three alleles present, any alleles beyond the top two that are above this threshold are reported.", optional = true)
    double minBasePct = 0.1;

    @Override
    public void onTraversalStart() {
        try {
            Utils.nonNull(outputFile);
            IOUtil.assertFileIsWritable(outputFile);
            outputStream = new PrintStream(outputFile);
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public void closeTool() {
        super.closeTool();

        outputStream.close();
    }

    protected static class SiteBaseCounter {
        public SimpleInterval interval;
        final List<FlaggedSite> flaggedSites;

        public SiteBaseCounter(SimpleInterval interval) {
            this.interval = interval;
            this.flaggedSites = new ArrayList<>();
        }

        public void merge(SiteBaseCounter counter2) {
            if (!interval.equals(counter2.interval)) {
                throw new RuntimeException(String.format("SiteBaseCounters do not have the same locus: %s: %d, %s: %d", interval.getContig(), interval.getStart(), counter2.interval.getContig(), counter2.interval.getStart()));
            }

            flaggedSites.addAll(counter2.flaggedSites);
        }

        public void doPrint(PrintStream out) {
            if (!flaggedSites.isEmpty()) {
                Set<String> totalSubjects = new HashSet<>();
                List<String> comments = new ArrayList<>();
                for (FlaggedSite f : flaggedSites) {
                    totalSubjects.add(f.sample);
                    comments.add(String.format("%s: [%s/%f/%d/%d]", f.sample, (char)(f.base.byteValue()), f.pct, f.baseCount, f.depth));
                }

                out.println(String.format("%s\t%d\t%d\t%s\t%s\t%s\t%s", interval.getContig(), interval.getStart()-1, interval.getEnd(), "MultiAllelicSite", totalSubjects.size(), "+", StringUtils.join(comments, ";")));
            }
        }

        public void addSite(FlaggedSite site) {
            flaggedSites.add(site);
        }
    }

    private static class FlaggedSite {
        public String sample;
        public Byte base;
        public double pct;
        public long depth;
        public long baseCount;

        public FlaggedSite(String sample, Byte base, double pct, long depth, long baseCount) {
            this.sample = sample;
            this.base = base;
            this.pct = pct;
            this.depth = depth;
            this.baseCount = baseCount;
        }
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        SiteBaseCounter counter = new SiteBaseCounter(referenceContext.getInterval());

        for (SAMReadGroupRecord rg: getHeaderForReads().getReadGroups()) {
            ReadPileup pileup = alignmentContext.getBasePileup().getPileupForSample(rg.getSample(), getHeaderForReads());
            int[] baseCounts = pileup.getBaseCounts();
            if (baseCounts == null) {
                continue;
            }

            int total = sumOfArray(baseCounts);
            if (total < minDepth) {
                continue;
            }

            Map<Integer, Double> baseCountPct = getBaseCountPctMap(baseCounts, total);
            if (baseCountPct.size() <= 2) {
                //nothing to do
            }
            else {
                List<Map.Entry<Integer, Double>> vals = sortByValue(baseCountPct);

                //accept the top two bases
                vals.remove(0);
                vals.remove(0);

                //iterate remaining by pct
                for (Map.Entry<Integer, Double> entry : vals) {
                    if (entry.getValue() >= minBasePct) {
                        try {
                            counter.addSite(new FlaggedSite(rg.getSample(), BaseUtils.BASES[entry.getKey()], entry.getValue(), total, baseCounts[entry.getKey()]));
                        }
                        catch (IndexOutOfBoundsException e) {
                            logger.warn(String.format("index of out bounds, %d, %d, %s, %d", entry.getKey(), baseCounts.length, referenceContext.getInterval().getContig(), referenceContext.getInterval().getStart()));
                            throw e;
                        }
                    }
                }
            }
        }

        counter.doPrint(outputStream);
    }

    public static <K, V extends Comparable<? super V>> List<Map.Entry<K, V>> sortByValue( Map<K, V> map ) {
        List<Map.Entry<K, V>> list =  new LinkedList<Map.Entry<K, V>>( map.entrySet() );
        Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
            public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 )
            {
                return (o1.getValue()).compareTo( o2.getValue() );
            }
        } );

        Collections.reverse(list);

        return list;
    }

    private Map<Integer, Double> getBaseCountPctMap(int[] baseCounts, int total)
    {
        Map<Integer, Double> ret = new HashMap<>();

        for (int i = 0; i < baseCounts.length; i++)
        {
            ret.put(i,  total == 0 ? 0 : ((double)baseCounts[i]) / total);
        }

        return ret;
    }

    private int sumOfArray(int[] array)
    {
        int ret = 0;
        for (Integer i : array)
        {
            ret += i;
        }

        return ret;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> ret = super.getDefaultReadFilters();
        ret.add(new ReadFilterLibrary.NotDuplicateReadFilter());
        return ret;
    }
}