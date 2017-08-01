package com.github.discvrseq.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.NumberFormat;
import java.util.*;

/**
 * Created by bimber on 7/31/2017.
 */
@CommandLineProgramProperties(oneLineSummary = "Provides genotyping summary for complex multi-genic loci, like KIR or MHC.", summary = ImmunoGenotyper.SUMMARY, programGroup = DiscvrSeqProgramGroup.class)
public class ImmunoGenotyper extends ReadWalker {
    protected static final String SUMMARY = "";

    @Argument(doc="File to which genotypes should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String out = null;

    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth.", optional = false)
    Integer minMappingQuality = 20;

    @Argument(fullName = "requireValidPair", shortName = "rvp", doc = "If true, only reads with a valid pair to the same reference will be considered", optional = true)
    boolean requireValidPair = false;


    @Argument(fullName = "minPctForRef", shortName = "minPctForRef", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, any reference with fewer than this fraction of reads mapping to it (of the total reads mapping) will be discarded.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minPctForRef = 0.01;

    @Argument(fullName = "minReadCountForRef", shortName = "minReadCountForRef", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, any reference with fewer than the provided number of reads mapping to it will be discarded.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minReadCountForRef = 0.01;

    @Argument(fullName = "minPctForLineageFiltering", shortName = "minPctForLineageFiltering", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, references will be grouped by lineage/allotype (based on the supplied file).  Within each group any hit sets below this threshold will be ignored.  Of the remaining, if any alleles are present across all groups, only those hits will be retained.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minPctForLineageFiltering = 0.01;

    @Argument(fullName = "referenceToLineageFile", shortName = "referenceToLineageFile", doc = "This is a simple tab-delimited file, no header, with two columns.  The first is the reference name, identical to that provided in the FASTA/BAM.  The second column is the lineage/allotype of this sequence.  This is used for grouping purposes and filtering.", optional = true)
    File referenceToLineageFile;

    @Argument(fullName = "minPctForExport", shortName = "minPctForExport", doc = "If provided, any genotype representing fewer than this fraction of mapped reads will be discarded.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minPctForExport = 0.01;

    @Argument(fullName = "minReadCountForExport", shortName = "minReadCountForExport", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, any genotype with fewer than the provided number of reads will be discarded.", optional = true, minValue = 0)
    Integer minReadCountForExport = 5;

    private ReferenceMatchTracker refTracker;
    private Map<String, String> nameToLineageMap;

    @Override
    public void onTraversalStart() {
        SAMFileHeader.SortOrder so = getHeaderForReads().getSortOrder();
        if (so != SAMFileHeader.SortOrder.queryname){
            throw new IllegalArgumentException("BAM must be in queryName sort order");
        }

        refTracker = new ReferenceMatchTracker();

        nameToLineageMap = new HashMap<>();
        if (referenceToLineageFile != null){
            try(BufferedReader reader = IOUtil.openFileForBufferedUtf8Reading(referenceToLineageFile)) {
                String line;
                while ((line = reader.readLine()) != null){
                    line = StringUtils.trimToNull(line);
                    if (line == null){
                        continue;
                    }

                    String[] cells = line.split("\t");
                    if (cells.length < 2){
                        throw new GATKException("Reference to Lineage/Allotype file must have two values per line");
                    }

                    if (nameToLineageMap.containsKey(cells[0])){
                        throw new GATKException("Each reference can only be provided once in the Lineage/Allotype file: " + cells[0]);
                    }

                    nameToLineageMap.put(cells[0], cells[1]);
                }
            }
            catch (IOException e){
                throw new GATKException(e.getMessage(), e);
            }
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new WellformedReadFilter());
    }

    private ReadAlignmentsTracker activeRead = null;

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (activeRead == null){
            activeRead = new ReadAlignmentsTracker(read.getName());
        }
        else if (!activeRead.activeReadName.equals(read.getName())) {
            refTracker.addRead(activeRead);

            activeRead = new ReadAlignmentsTracker(read.getName());
        }

        if (read.getMappingQuality() < minMappingQuality){
            activeRead.addLowMapqAlignment();
            return;
        }

        //determine if has mismatches
        if (!read.hasAttribute("NM")) {
            throw new IllegalArgumentException("Read lacks NM tag: " + read.getName());
        }

        if (read.getAttributeAsInteger("NM") > 0){
            activeRead.addMismatchAlignment();
            return;
        }

        activeRead.addAlignment(read);

        return;
    }

    @Override
    public Object onTraversalSuccess() {
        //outputWriter = new PrintWriter(OUTPUT);
        if (activeRead != null){
            refTracker.addRead(activeRead);
        }

        //perform filtering of hits
        filterByReference(refTracker);
        filterByLineage(refTracker);

        //write out tables
        Path outputFile = IOUtils.getPath(out);
        NumberFormat numberFormat = NumberFormat.getNumberInstance();
        numberFormat.setMinimumFractionDigits(3);

        IOUtil.assertFilesAreWritable(Arrays.asList(outputFile.toFile()));
        try (PrintWriter outWriter = new PrintWriter(IOUtil.openFileForBufferedWriting(outputFile.toFile()))){
            outWriter.println(StringUtils.join(Arrays.asList("RefNames", "Lineage/Allotypes", "TotalReads", "PercentOfTotal", "PercentOfTotalIncludingUnmapped"), "\t"));

            for (HitSet hs : refTracker.hitMap.values()) {
                Double pct = hs.readNames.size() / (double) refTracker.readPairsWithHits;
                Double pct2 = hs.readNames.size() / (double)(refTracker.readPairsWithHits + refTracker.readPairsNoHits);

                if (hs.readNames.size() < minReadCountForExport){
                    continue;
                }

                if (pct < minPctForExport){
                    continue;
                }

                outWriter.println(StringUtils.join(Arrays.asList(
                        StringUtils.join(hs.refNames, ","),
                        StringUtils.join(hs.getLineages(nameToLineageMap), ","),
                        String.valueOf(hs.readNames.size()),
                        numberFormat.format(pct),
                        numberFormat.format(pct2)
                ), "\t"));
            }
        }

        return super.onTraversalSuccess();
    }

    private class ReferenceMatchTracker {
        private Map<String, HitSet> hitMap = new HashMap<>();
        private int readPairsWithHits = 0;
        private int readPairsNoHits = 0;
        private int totalReadsFailedForMapq = 0;
        private int totalReadsFailedForValidPair = 0;

        public ReferenceMatchTracker(){

        }

        public void addRead(ReadAlignmentsTracker tracker){
            Set<String> hits = tracker.getHits(requireValidPair);
            if (!hits.isEmpty()){
                readPairsWithHits++;

                String key = HitSet.getKey(hits);
                if (!hitMap.containsKey(key)){
                    hitMap.put(key, new HitSet(hits));
                }

                hitMap.get(key).addRead(tracker);
            }
            else {
                readPairsNoHits++;

                if (tracker.lowMapqAlignments > 0){
                    totalReadsFailedForMapq++;
                }

                if (requireValidPair && (!tracker.forwardPerfectHits.isEmpty() || !tracker.reversePerfectHits.isEmpty())){
                    totalReadsFailedForValidPair++;
                }
            }
        }
    }

    private class ReadAlignmentsTracker {
        private String activeReadName;

        private Set<String> forwardPerfectHits = new HashSet<>();
        private Set<String> reversePerfectHits = new HashSet<>();
        private int lowMapqAlignments = 0;
        private int mismatchAlignments = 0;

        public ReadAlignmentsTracker(String activeReadName){
            this.activeReadName = activeReadName;
        }

        public void addAlignment(GATKRead record) {
            if (!record.isPaired() || record.isFirstOfPair()){
                forwardPerfectHits.add(record.getContig());
            }
            else if (record.isSecondOfPair()){
                reversePerfectHits.add(record.getContig());
            }
        }

        public void addLowMapqAlignment() {
            lowMapqAlignments++;
        }

        public void addMismatchAlignment() {
            mismatchAlignments++;
        }

        public Set<String> getHits(boolean requireValidPair)
        {
            Set<String> ret = new HashSet<>();
            ret.addAll(forwardPerfectHits);
            if (requireValidPair){
                ret.retainAll(reversePerfectHits);
            }
            else {
                ret.addAll(reversePerfectHits);
            }

            return ret;
        }
    }

    private static class HitSet {
        private Set<String> readNames = new HashSet<>();
        private Set<String> refNames = new TreeSet<>();

        private int forward = 0;
        private int reverse = 0;
        private int valid_pair = 0;

        public HitSet(Collection<String> refNames) {
            this.refNames.addAll(refNames);
        }

        public void append(HitSet other) {
            forward += other.forward;
            reverse += other.reverse;
            valid_pair += other.valid_pair;
            readNames.addAll(other.readNames);
        }

        public String getKey() {
            return getKey(refNames);
        }

        public void addRead(ReadAlignmentsTracker tracker){
            this.readNames.add(tracker.activeReadName);
        }

        public static String getKey(Collection<String> refNames) {
            List<String> ret = new ArrayList<>(refNames);
            Collections.sort(ret);

            return StringUtils.join(ret, "||");
        }

        public Set<String> getLineages(Map<String, String> referenceToLineageMap){
            TreeSet ret = new TreeSet();
            for (String refName : refNames){
                ret.add(referenceToLineageMap.containsKey(refName) ? referenceToLineageMap.get(refName) : refName);
            }

            return ret;
        }
    }

    private void filterByReference(ReferenceMatchTracker refTracker){
        int readsHelpedByAlleleFilters = 0;

        //build total by ref
        Map<String, Integer> totalByReference = new HashMap<>();
        int totalReads = 0;
        for (HitSet hs : refTracker.hitMap.values()) {
            for (String refName : hs.refNames) {
                int total = totalByReference.containsKey(refName) ? totalByReference.get(refName) : 0;
                total += hs.readNames.size();
                totalByReference.put(refName, total);
            }

            totalReads += hs.readNames.size();
        }

        //make blacklist
        Set<String> disallowedReferences = new HashSet<>();
        for (String refName : totalByReference.keySet()){
            int totalForRef = totalByReference.get(refName);
            double pct = ((double) totalByReference.get(refName) / totalReads);

            if (minReadCountForRef != null && totalForRef < minReadCountForRef) {
                //_skippedReferencesByRead++;
                //getLogger().debug("Reference discarded due to read count: " + refName + " / " + distinctStageTwoReads + " / " + totalForRef + " / " + pct + "%");
                //msg = "**skipped due to read count";
                disallowedReferences.add(refName);
            }
            else if (minPctForRef != null && pct < minPctForRef) {
                //_skippedReferencesByPct++;
                //getLogger().debug("Reference discarded due to percent: " + refName + " / " + distinctStageTwoReads + " / " + totalForRef + " / " + pct + "%");
                //msg = "**skipped due to percent";
                disallowedReferences.add(refName);
            }
        }

        //now actually filter:
        //then actually use these for filtering
        Map<String, HitSet> newHitMap = new HashMap<>();
        for (HitSet hs : refTracker.hitMap.values()) {
            List<String> refNames = new ArrayList<>(hs.refNames);
            refNames.removeAll(disallowedReferences);
            if (refNames.isEmpty()) {
                //refTracker.unaligned.addAll(hs.readNames);
            }
            else {
                if (refNames.size() != hs.refNames.size()) {
                    readsHelpedByAlleleFilters++;
                }

                //merge sets
                String newKey = HitSet.getKey(refNames);
                HitSet hs2 = newHitMap.containsKey(newKey) ? newHitMap.get(newKey) : new HitSet(refNames);
                hs2.append(hs);
                newHitMap.put(newKey, hs2);
            }
        }

        refTracker.hitMap = newHitMap;
    }

    private void filterByLineage(ReferenceMatchTracker refTracker){
        if (nameToLineageMap.isEmpty()){
            logger.info("no reference to lineage/allotype file provided, cannot perform filtering");
            return;
        }

        Map<String, HitSet> newHitSetMap = new HashMap<>();

        //build a map of distinct sets by lineage
        Map<String, List<HitSet>> resultByLineage = new HashMap<>();
        Map<String, Integer> totalByLineage = new HashMap<>();
        for (String key : refTracker.hitMap.keySet()) {
            HitSet hs = refTracker.hitMap.get(key);
            Set<String> distinctLineages = new HashSet<>();
            for (String refName: hs.refNames) {
                if (!nameToLineageMap.containsKey(refName)) {
                    //if we have missing lineages, abort and keep data as-is
                    distinctLineages.clear();
                    break;
                }

                distinctLineages.add(nameToLineageMap.get(refName));
            }

            if (distinctLineages.size() == 1) {
                String lineage = distinctLineages.iterator().next();
                if (!resultByLineage.containsKey(lineage)) {
                    resultByLineage.put(lineage, new ArrayList<>());
                    totalByLineage.put(lineage, 0);
                }

                resultByLineage.get(lineage).add(hs);
                totalByLineage.put(lineage, totalByLineage.get(lineage) + hs.readNames.size());
            }
            else {
                newHitSetMap.put(key, hs);
            }
        }

        //now filter by lineage
        logger.info("total lineages being inspected: " + resultByLineage.size());
        for (String lineage : resultByLineage.keySet()) {
            List<HitSet> sets = resultByLineage.get(lineage);
            if (sets.size() == 1) {
                newHitSetMap.put(sets.get(0).getKey(), sets.get(0));
                continue;
            }

            if (!totalByLineage.containsKey(lineage)) {
                logger.error("unable to find lineage, cannot filter: [" + lineage + "]");
                for (HitSet hs : sets) {
                    newHitSetMap.put(hs.getKey(), hs);
                }
                continue;
            }

            Set<String> sharedRefNames = null;
            int setsSkipped = 0;
            for (HitSet hs : sets) {
                double pctOfLineage = (double)hs.readNames.size() / (double)totalByLineage.get(lineage);
                if (pctOfLineage < minPctForLineageFiltering) {
                    setsSkipped++;
                    continue;
                }

                if (sharedRefNames == null) {
                    sharedRefNames = new HashSet<>();
                    sharedRefNames.addAll(hs.refNames);
                }
                else {
                    sharedRefNames.retainAll(hs.refNames);
                    if (sharedRefNames.isEmpty()){
                        break;
                    }
                }
            }

            logger.debug("total sets skipped due to pct: " + setsSkipped);
            if (sharedRefNames == null || sharedRefNames.isEmpty()) {
                //if empty, there are no alleles common to all, so keep original data
                for (HitSet hs : sets) {
                    newHitSetMap.put(hs.getKey(), hs);
                }
            }
            else {
                //merge and make new
                HitSet merged = new HitSet(sharedRefNames);
                for (HitSet hs : sets) {
                    //if below the threshold, leave as is
                    double pctOfLineage = (double)hs.readNames.size() / (double)totalByLineage.get(lineage);
                    if (pctOfLineage < minPctForLineageFiltering) {
                        if (newHitSetMap.containsKey(hs.getKey())) {
                            newHitSetMap.get(hs.getKey()).append(hs);
                        }
                        else {
                            newHitSetMap.put(hs.getKey(), hs);
                        }

                        continue;
                    }

                    merged.append(hs);
                }

                if (newHitSetMap.containsKey(merged.getKey())) {
                    newHitSetMap.get(merged.getKey()).append(merged);
                }
                else {
                    newHitSetMap.put(merged.getKey(), merged);
                }
            }
        }

        refTracker.hitMap = newHitSetMap;
    }
}
