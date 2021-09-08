package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqDevProgramGroup;
import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
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
 * This tool will generate genotype calls for complex loci (such as many immune genes), from next-generation sequence data.
 * In general, most next-generation sequencing tasks (WGS, RNA-seq, etc.) expect to align each read uniquely in the genome and frequently discard or ignore
 * those reads that do not. While this is sufficient for much of the genome, complex, polygenic regions are frequently systematically filtered or missed by many
 * variant-calling or genotyping methods. Further, the reference genome often has a poor representation of complex loci, including those that are polygenic or have copy-number differences
 * in the population (often the reference genome has a single copy of the gene), or situations where the set of genes varies widely between individuals (in which case a single reference genome
 * poorly captures variation). This tool applies the same general approach as originally published for macaque MHC genotyping (<a href="https://www.ncbi.nlm.nih.gov/pubmed/19820716">https://www.ncbi.nlm.nih.gov/pubmed/19820716</a>),
 * though is a much newer implementation. This tool can either be run using targeted sequencing (i.e. amplicon), or using whole genome DNA- or RNA-seq data.
 *
 * To use this tool, first align your raw reads to a specialized database. For MHC, KIR, or NCR genotyping, we build a specialized database of alleles (<a href="https://www.ebi.ac.uk/ipd/">IPD is a good starting point</a>).
 * We often use an aligner that will retain multiple matches per read (<a href="https://github.com/wanpinglee/MOSAIK">MOSAIK</a>, for example); however, we use BWA-mem in some cases as well. This BAM is the input
 * for the this tool. ImmunoGenotyper iterates this BAM, capturing the set of contigs each read aligned against, within the provided number of mismatches. We expect ambiguity in many cases, with reads aligning
 * to groups of related alleles (i.e. A002*001:02 and A*002*001:01). Finally, we series of filters can be run (described for each argument below) that attempt to collapse these genotypes. A final table
 * is produced. Note: the reference data is essential for this tool to operate as-expected, and thoughtful design of that database for your application is critical.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar ImmunoGenotyper \
 *     -R reference.fasta \
 *     -I myBam.bam \
 *     -mm 1 \
 *     -minReadCountForRef 10 \
 *     -O outputBaseDir
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Provides genotyping summary for complex multigenic loci, like KIR or MHC.",
        summary = ImmunoGenotyper.SUMMARY,
        programGroup = DiscvrSeqDevProgramGroup.class)
public class ImmunoGenotyper extends ReadWalker {
    protected static final String SUMMARY = "Genotyping tool for complex loci";

    public static final String GENOTYPE_EXTENSION = ".genotypes.txt";
    public static final String MISMATCH_EXTENSION = ".mismatches.txt";
    public static final String SUMMARY_EXTENSION = ".summary.txt";

    @Argument(doc="Prefix for output files", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outPrefix = null;

    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth.", optional = false)
    Integer minMappingQuality = 0;

    @Argument(fullName = "requireValidPair", shortName = "rvp", doc = "If true, only reads with a valid pair to the same reference will be considered", optional = true)
    boolean requireValidPair = false;

    @Argument(fullName = "minPctForRef", shortName = "minPctForRef", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, any reference with fewer than this fraction of reads mapping to it (of the total reads mapping) will be discarded.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minPctForRef = 0.01;

    @Argument(fullName = "minReadCountForRef", shortName = "minReadCountForRef", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, any reference with fewer than the provided number of reads mapping to it will be discarded.", optional = true, minValue = 0)
    Integer minReadCountForRef = 5;

    @Argument(fullName = "minPctForLineageFiltering", shortName = "minPctForLineageFiltering", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, references will be grouped by lineage/allotype (based on the supplied file).  Within each group any hit sets below this threshold will be ignored.  Of the remaining, if any alleles are present across all groups, only those hits will be retained.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minPctForLineageFiltering = 0.01;

    @Argument(fullName = "referenceToLineageFile", shortName = "referenceToLineageFile", doc = "This is a simple tab-delimited file, no header, with two columns.  The first is the reference name, identical to that provided in the FASTA/BAM.  The second column is the lineage/allotype of this sequence.  This is used for grouping purposes and filtering.", optional = true)
    File referenceToLineageFile;

    @Argument(fullName = "minPctForExport", shortName = "minPctForExport", doc = "If provided, any genotype representing fewer than this fraction of mapped reads will be discarded.", optional = true, maxValue = 1.0, minValue = 0.0)
    Double minPctForExport = 0.00;

    @Argument(fullName = "minReadCountForExport", shortName = "minReadCountForExport", doc = "This is part of the filtering strategy for ambiguous hits.  If provided, any genotype with fewer than the provided number of reads will be discarded.", optional = true, minValue = 0)
    Integer minReadCountForExport = 5;

    @Argument(fullName = "mismatchesTolerated", shortName = "mm", doc = "The maximum number of mismatches tolerated for alignment.  If a given read has multiple alignments, those with the fewest mismatches will be kept (irrespective of length).", optional = true, minValue = 0)
    Integer mismatchesTolerated = 0;

    @Argument(fullName = "minAlignmentLength", shortName = "minAlignmentLength", doc = "Alignments shorted than this value will be discarded.", optional = true, minValue = 0)
    Integer minAlignmentLength = 40;

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

    private class Filter extends ReadFilter {
        private static final long serialVersionUID = 1L;

        private ReadFilter filter = null;

        public Filter() {
        }

        @Override
        public void setHeader(SAMFileHeader header) {
            super.setHeader(header);
            createFilter();
        }

        public Filter(final SAMFileHeader header) {
            setHeader(header);
        }

        private void createFilter() {
            final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(samHeader);

            filter = ReadFilterLibrary.VALID_ALIGNMENT_START
                    .and(ReadFilterLibrary.VALID_ALIGNMENT_END);
                    //.and(alignmentAgreesWithHeader)
                    //.and(ReadFilterLibrary.HAS_READ_GROUP)
                    //.and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                    //.and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                    //.and(ReadFilterLibrary.SEQ_IS_STORED)
                    //.and(ReadFilterLibrary.CIGAR_CONTAINS_NO_N_OPERATOR);
        }

        @Override
        public boolean test(GATKRead read) {
            return filter.test(read);
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new Filter());
    }

    private ReadAlignmentsTracker activeRead = null;

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (read.isUnmapped()){
            return;
        }

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

        activeRead.addAlignment(read, read.getAttributeAsInteger("NM"));
    }

    @Override
    public Object onTraversalSuccess() {
        //finalize the last group
        if (activeRead != null){
            refTracker.addRead(activeRead);
        }

        NumberFormat numberFormat = NumberFormat.getNumberInstance();
        numberFormat.setMinimumFractionDigits(3);

        //perform filtering of hits
        List<String> messages = new ArrayList<>();
        double total = (double)refTracker.readPairsNoHits + refTracker.readPairsWithHits;
        messages.add("Read pairs with hits: " + refTracker.readPairsWithHits + " (" + numberFormat.format(refTracker.readPairsWithHits / total) + ")");
        messages.add("Read pairs without hits: " + refTracker.readPairsNoHits + " (" + numberFormat.format(refTracker.readPairsNoHits / total) + ")");
        messages.add("Failed due to MAPQ: " + refTracker.totalReadsFailedForMapq + " (" + numberFormat.format(refTracker.totalReadsFailedForMapq / total) + ")");
        messages.add("Failed due to length: " + refTracker.totalReadsFailedForLength + " (" + numberFormat.format(refTracker.totalReadsFailedForLength / total) + ")");
        messages.add("Failed due to no valid pair: " + refTracker.totalReadsFailedForValidPair + " (" + numberFormat.format(refTracker.totalReadsFailedForValidPair / total) + ")");
        messages.add("Failed due to mismatches: " + refTracker.totalAlignmentsFailedForMismatch + " (" + numberFormat.format(refTracker.totalAlignmentsFailedForMismatch / total) + ")");

        filterByReference(refTracker, messages);
        filterByLineage(refTracker, messages);

        //write out tables
        Path outputFile = IOUtils.getPath(outPrefix + GENOTYPE_EXTENSION);
        IOUtil.assertFilesAreWritable(Arrays.asList(outputFile.toFile()));
        try (PrintWriter outWriter = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(outputFile.toFile()))){
            outWriter.println(StringUtils.join(Arrays.asList("RefNames", "Lineage/Allotypes", "TotalReads", "PercentOfTotal", "PercentOfTotalIncludingUnmapped"), "\t"));

            messages.add("Exporting final groups:");

            int groupsSkipped = 0;
            for (String key : refTracker.hitMap.keySet()) {
                HitSet hs = refTracker.hitMap.get(key);
                Double pct = hs.readNames.size() / (double) refTracker.readPairsWithHits;
                Double pct2 = hs.readNames.size() / (double)(refTracker.readPairsWithHits + refTracker.readPairsNoHits);

                if (hs.readNames.size() < minReadCountForExport){
                    messages.add("Discarded due to count: " + key + " / " + hs.readNames.size() + " / " + numberFormat.format(pct));
                    groupsSkipped++;
                    continue;
                }

                if (pct < minPctForExport){
                    messages.add("Discarded due to percent: " + key + " / " + hs.readNames.size() + " / " + numberFormat.format(pct));
                    groupsSkipped++;
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

            messages.add("Groups skipped due to low read count or percent: " + groupsSkipped);
        }

        Path summaryFile = IOUtils.getPath(outPrefix + SUMMARY_EXTENSION);
        IOUtil.assertFilesAreWritable(Arrays.asList(summaryFile.toFile()));
        try (PrintWriter outWriter = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(summaryFile.toFile()))){
            messages.forEach(message -> outWriter.println(message));
        }

        //Also log to console:
        messages.forEach(message -> logger.info(message));

        Path mismatchFile = IOUtils.getPath(outPrefix + MISMATCH_EXTENSION);
        IOUtil.assertFilesAreWritable(Arrays.asList(mismatchFile.toFile()));
        try (PrintWriter outWriter = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(mismatchFile.toFile()))){
            outWriter.println(StringUtils.join(Arrays.asList("RefName", "TotalReads", "ReasonForFailure"), "\t"));
            for (String refName : refTracker.mismatchMap.keySet()){
                AlignmentMismatch am = refTracker.mismatchMap.get(refName);
                outWriter.println(StringUtils.join(Arrays.asList(
                    refName,
                    String.valueOf(am.totalReads),
                    StringUtils.join(am.reasonsForFailure, ",")
                ), "\t"));
            }
        }

        return super.onTraversalSuccess();
    }

    private class AlignmentMismatch {
        String refName;
        int totalReads = 0;
        Set<String> reasonsForFailure = new HashSet<>();

        AlignmentMismatch(String refName){
            this.refName = refName;
        }

        private void addRead(String reason){
            totalReads++;
            reasonsForFailure.add(reason);
        }
    }

    private class ReferenceMatchTracker {
        private Map<String, HitSet> hitMap = new HashMap<>();
        private Map<String, AlignmentMismatch> mismatchMap = new HashMap<>();
        private int readPairsWithHits = 0;
        private int readPairsNoHits = 0;
        private int totalReadsFailedForMapq = 0;
        private int totalReadsFailedForValidPair = 0;
        private int totalReadsFailedForLength = 0;
        private int totalAlignmentsFailedForMismatch = 0;

        public ReferenceMatchTracker(){

        }

        public void addRead(ReadAlignmentsTracker tracker){
            ReadHit hit = tracker.getHits(requireValidPair);
            if (!hit.hits.isEmpty()){
                readPairsWithHits++;

                String key = HitSet.getKey(hit.hits);
                if (!hitMap.containsKey(key)){
                    hitMap.put(key, new HitSet(hit.hits));
                }

                hitMap.get(key).addRead(tracker, hit.hasForward, hit.hasReverse);
            }
            else {
                readPairsNoHits++;

                if (tracker.lowMapqAlignments > 0){
                    totalReadsFailedForMapq++;
                }

                if (tracker.shortAlignments.size() > 0){
                    totalReadsFailedForLength++;
                }

                if (tracker.mismatchAlignments.size() > 0){
                    totalAlignmentsFailedForMismatch++;
                    for (String refName : tracker.mismatchAlignments){
                        if (!mismatchMap.containsKey(refName)){
                            mismatchMap.put(refName, new AlignmentMismatch(refName));
                        }

                        mismatchMap.get(refName).addRead("Mismatches");
                    }
                }

                if (requireValidPair && (!tracker.forwardPerfectHits.isEmpty() || !tracker.reversePerfectHits.isEmpty())){
                    totalReadsFailedForValidPair++;

                    ReadHit hitNoValidPair = tracker.getHits(false);
                    if (!hitNoValidPair.hits.isEmpty()){
                        for (String refName : hitNoValidPair.hits){
                            if (!mismatchMap.containsKey(refName)){
                                mismatchMap.put(refName, new AlignmentMismatch(refName));
                            }

                            mismatchMap.get(refName).addRead("NoValidPair");
                        }
                    }
                }
            }
        }
    }

    private class ReadHit {
        Set<String> hits = new HashSet<>();
        boolean hasForward;
        boolean hasReverse;
    }

    private class ReadAlignmentsTracker {
        private String activeReadName;

        private Map<Integer, Set<String>> forwardPerfectHits = new HashMap<>();
        private Map<Integer, Set<String>> reversePerfectHits = new HashMap<>();
        private int lowMapqAlignments = 0;
        private Set<String> mismatchAlignments = new HashSet<>();
        private Set<String> shortAlignments = new HashSet<>();

        public ReadAlignmentsTracker(String activeReadName){
            this.activeReadName = activeReadName;
        }

        public void addAlignment(GATKRead record, int nm) {
            if (record.getLength() < minAlignmentLength){
                shortAlignments.add(record.getContig());
            }
            else if (nm > mismatchesTolerated){
                mismatchAlignments.add(record.getContig());
            }
            else {
                if (!forwardPerfectHits.containsKey(nm)){
                    forwardPerfectHits.put(nm, new HashSet<>());
                    reversePerfectHits.put(nm, new HashSet<>());
                }

                if (!record.isPaired() || record.isFirstOfPair()){
                    forwardPerfectHits.get(nm).add(record.getContig());
                }
                else if (record.isSecondOfPair()) {
                    reversePerfectHits.get(nm).add(record.getContig());
                }
            }
        }

        public void addLowMapqAlignment() {
            lowMapqAlignments++;
        }

        public ReadHit getHits(boolean requireValidPair)
        {
            if (forwardPerfectHits.isEmpty()){
                return new ReadHit();
            }

            Integer lowestForward = new TreeSet<>(forwardPerfectHits.keySet()).iterator().next();
            Integer lowestReverse = new TreeSet<>(reversePerfectHits.keySet()).iterator().next();
            Set<String> hits = new HashSet<>();
            hits.addAll(forwardPerfectHits.get(lowestForward));
            if (requireValidPair){
                hits.retainAll(reversePerfectHits.get(lowestReverse));
                ReadHit ret = new ReadHit();
                ret.hits.addAll(hits);
                ret.hasForward = true;
                ret.hasReverse = true;
                return ret;

            }
            else {
                hits.addAll(reversePerfectHits.get(lowestReverse));
                ReadHit ret = new ReadHit();
                ret.hits.addAll(hits);
                ret.hasForward = !forwardPerfectHits.get(lowestForward).isEmpty();
                ret.hasReverse = !reversePerfectHits.get(lowestReverse).isEmpty();
                return ret;
            }
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

        public void addRead(ReadAlignmentsTracker tracker, boolean isForward, boolean isReverse){
            this.readNames.add(tracker.activeReadName);

            if (isForward){
                forward++;
            }

            if (isReverse){
                reverse++;
            }
        }

        public static String getKey(Collection<String> refNames) {
            List<String> ret = new ArrayList<>(refNames);
            Collections.sort(ret);

            return StringUtils.join(ret, "||");
        }

        public Set<String> getLineages(Map<String, String> referenceToLineageMap){
            TreeSet<String> ret = new TreeSet<>();
            for (String refName : refNames){
                ret.add(referenceToLineageMap.getOrDefault(refName, refName));
            }

            return ret;
        }
    }

    private void filterByReference(ReferenceMatchTracker refTracker, List<String> messages){
        messages.add("Filtering by reference:");
        int readsHelpedByAlleleFilters = 0;

        //build total by ref
        Map<String, Integer> totalByReference = new HashMap<>();
        int totalReads = 0;
        for (HitSet hs : refTracker.hitMap.values()) {
            for (String refName : hs.refNames) {
                int total = totalByReference.getOrDefault(refName, 0);
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
                messages.add("Discarded due to read count: " + refName + " / " + totalForRef + " / " + pct);
                disallowedReferences.add(refName);
            }
            else if (minPctForRef != null && pct < minPctForRef) {
                messages.add("Discarded due to percent: " + refName + " / " + totalForRef + " / " + pct);
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
                //messages.add("Set condensed: " + StringUtils.join(hs.refNames, ",") + " -> None");
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

        messages.add("Groups before/after filtering by reference: " + refTracker.hitMap.size() + "/" + newHitMap.size());

        refTracker.hitMap = newHitMap;
    }

    private void filterByLineage(ReferenceMatchTracker refTracker, List<String> messages){
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
                    sharedRefNames = new HashSet<>(hs.refNames);
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

        messages.add("Groups before/after filtering by lineage/allotype: " + refTracker.hitMap.size() + "/" + newHitSetMap.size());

        refTracker.hitMap = newHitSetMap;
    }
}
