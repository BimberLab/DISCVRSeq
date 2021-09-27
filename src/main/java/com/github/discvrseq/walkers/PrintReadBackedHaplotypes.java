package com.github.discvrseq.walkers;


import com.github.discvrseq.tools.DiscvrSeqDevProgramGroup;
import com.github.discvrseq.util.CigarPositionIterable;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.IntervalWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * This tool will extract the reads from a BAM over each of the provided intervals, reconstruct the local haplotypes (using simple logic and only relying on regions with coverage), and
 * produce a table listing the frequency of every unique haplotype.  It was originally created to inspect amplicon-based deep sequencing, such as evaluating CRISPR edits.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadBackedHaplotypes \
 *     -R reference.fasta \
 *     -L chr01:100-200 \
 *     -I myBam.bam \
 *     -O output.txt
 * </pre>
 * <h3>Usage example, <a href="https://gatkforums.broadinstitute.org/gatk/discussion/1319/collected-faqs-about-interval-lists"></a>supplying intervals from a file:</a></h3>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadBackedHaplotypes \
 *     -R reference.fasta \
 *     -L myFile.intervals \
 *     -I myBam.bam \
 *     -O output.txt
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Reconstructs and prints a list of distinct read-backed sequences over the supplied intervals.",
        oneLineSummary = "Reconstructs and prints a list of distinct read-backed sequences over the supplied intervals",
        programGroup = DiscvrSeqDevProgramGroup.class
)
public class PrintReadBackedHaplotypes extends IntervalWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    @Argument(fullName = "requiredCoverageFraction", shortName = "rc", doc = "If provided, only reads or read pairs with complete coverage over at least this fraction of the interval will be included")
    private Double requiredCoverageFraction  = 0.0;

    @Argument(fullName = "minQual", shortName = "mq", doc = "If specified, bases with quality lower than this value will be converted to N")
    final int minQual = 0;

    @Argument(fullName = "minMappingQual", shortName = "mmq", doc = "If specified, only alignments with mapping quality above this value will be included")
    final int minMappingQual = 10;

    private PrintStream outputStream = null;

    private SamReader bamReader = null;
    private ReadFilter readFilter;

    @Override
    public void onTraversalStart() {
        if (!hasUserSuppliedIntervals()) {
            throw new UserException.BadInput("Must supply a list of intervals on the command line");
        }

        if (readArguments.getReadPaths().size() > 1) {
            throw new UserException.BadInput("Only one BAM at a time is currently supported");
        }

        try {
            Utils.nonNull(outputFile);
            IOUtil.assertFileIsWritable(outputFile);
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }

        SamReaderFactory fact = SamReaderFactory.makeDefault();
        Path bam = readArguments.getReadPaths().get(0);
        bamReader = fact.open(bam);

        readFilter = ReadFilter.fromList(getDefaultReadFilters(), getHeaderForReads());
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private long failedFilters = 0L;
    private long readsInspected = 0L;
    private long totalDroppedForCoverage = 0L;
    private Map<Integer, Long> readTotalHist = new TreeMap<>();

    private final Map<SimpleInterval, Map<Character[][], Integer>> resultMap = new HashMap<>();

    @Override
    public void apply(SimpleInterval interval, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        Map<Character[][], Integer> results = new HashMap<>();
        Map<String, List<SAMRecord>> readMap = queryOverlappingReads(interval);
        for (String readName : readMap.keySet()) {
            readsInspected++;

            Character[][] arr = processGroup(interval, readMap.get(readName));
            if (arr == null) {
                continue;
            }

            boolean found = false;
            for (Character[][] val : results.keySet())
            {
                if (Arrays.deepEquals(val, arr))
                {
                    results.put(val, results.get(val) + 1);
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                results.put(arr, 1);
            }
        }

        resultMap.put(interval, results);
    }

    private Map<String, List<SAMRecord>> queryOverlappingReads(SimpleInterval interval) {
        Map<String, List<SAMRecord>> readMap = new HashMap<>();
        try (SAMRecordIterator it = bamReader.queryOverlapping(interval.getContig(), interval.getStart(), interval.getEnd())) {
            while (it.hasNext()) {
                GATKRead read = new SAMRecordToGATKReadAdapter(it.next());
                if (!readFilter.test(read)){
                    failedFilters++;
                    continue;
                }

                List<SAMRecord> list = readMap.getOrDefault(read.getName(), new ArrayList<>());
                list.add(read.convertToSAMRecord(getHeaderForReads()));
                readMap.put(read.getName(), list);
            }
        }

        return readMap;
    }

    private static final int MAX_NON_COVER_WINDOW = 200;

    private Character[][] processGroup(SimpleInterval interval, List<SAMRecord> reads) {
        Character[][] arr = new Character[interval.size()][];
        Integer[][] qualArr = new Integer[interval.size()][];

        readTotalHist.put(reads.size(), readTotalHist.getOrDefault(reads.size(), 0L) + 1);

        reads.forEach(read -> processRead(read, interval, arr, qualArr));

        //NOTE: for a given read group, large deletions can appear as an internal region w/o coverage
        boolean encounteredCoverage = false;
        for (int idx = 0;idx < arr.length; idx++){
            if (arr[idx] != null) {
                encounteredCoverage = true;
                continue;
            }

            if (!encounteredCoverage) {
                continue;
            }

            int nextNotNummIndex = IntStream.range(idx + 1, Math.min(idx + MAX_NON_COVER_WINDOW, arr.length))
                            .filter(i -> arr[i] != null)
                            .findFirst().orElse(-1);

            if (nextNotNummIndex > 0) {
                arr[idx] = new Character[]{'-'};
            }
        }

        if (requiredCoverageFraction > 0) {
            long totalCovered = IntStream.range(0, arr.length)
                    .filter(i -> arr[i] != null)
                    .count();

            double fraction = ((double)totalCovered / arr.length);
            if (fraction < requiredCoverageFraction) {
                totalDroppedForCoverage++;
                return null;
            }
        }

        return arr;
    }

    private void processRead(SAMRecord r, SimpleInterval interval, Character[][] arr, Integer[][] qualArr)
    {
        //add this value to a reference coordinate to find array position
        final int offset = interval.getStart() * -1;

        CigarPositionIterable cpi = new CigarPositionIterable(r);
        CigarPositionIterable.CigarIterator ci = cpi.iterator();

        int effectiveInsertIdx = 0;
        while (ci.hasNext())
        {
            CigarPositionIterable.PositionInfo pi = ci.next();

            //note: getRefPosition is 0-based, interval is 1-based
            if (pi.getRefPosition() + 1 < interval.getStart())
            {
                continue;
            }
            else if (pi.isSkipped())
            {
                continue;
            }
            else if (pi.getRefPosition() + 1 > interval.getEnd())
            {
                break;
            }

            if (pi.getInsertIndex() == 0)
            {
                //reset each time we hit a new position
                effectiveInsertIdx = 0;
            }

            //getRefPosition() is 0-based
            int arrayPos = pi.getRefPosition() + offset + 1;


            if (pi.isIndel())
            {
                if (pi.isDel())
                {
                    if (arr[arrayPos] == null)
                    {
                        arr[arrayPos] = new Character[]{pi.getBaseQuality() < minQual ? 'N' : '-'};
                        qualArr[arrayPos] = new Integer[]{pi.getBaseQuality()};
                    }
                    else
                    {
                        mergePositions(arr, qualArr, arrayPos, 0, pi.getBaseQuality() < minQual ? 'N' : '-', pi.getBaseQuality(), pi);
                    }
                }
                else if (pi.isInsertion() && pi.getBaseQuality() >= minQual)
                {
                    logger.info("indel: " + pi.getRecord().getReadName() + ", " + pi.getRefPosition());
                    //TODO: account for second mate
                    effectiveInsertIdx++;
                    Character[] posArr = arr[arrayPos];
                    if (posArr == null)
                    {
                        throw new IllegalArgumentException("No previous array for position: " + pi.getRefPosition());
                    }

                    Character[] newArr = Arrays.copyOf(posArr, effectiveInsertIdx + 1);
                    newArr[effectiveInsertIdx] = (char)pi.getReadBase();

                    arr[arrayPos] = newArr;
                }
            }
            else
            {
                if (arr[arrayPos] == null)
                {
                    arr[arrayPos] = new Character[]{pi.getBaseQuality() < minQual ? 'N' : (char) pi.getReadBase()};
                    qualArr[arrayPos] = new Integer[]{pi.getBaseQuality()};
                }
                else
                {
                    mergePositions(arr, qualArr, arrayPos, 0, pi.getBaseQuality() < minQual ? 'N' : (char) pi.getReadBase(), pi.getBaseQuality(), pi);
                }
            }
        }
    }

    private void mergePositions(Character[][] arr, Integer[][] qualArray, int arrayPos, int idx, char base, int qual, CigarPositionIterable.PositionInfo pi)
    {
        char existing = Character.toUpperCase(arr[arrayPos][idx]);
        if (existing == 'N')
        {
            arr[arrayPos][idx] = base;
        }
        else if (base == 'N')
        {
            return;
        }
        else if (existing != base)
        {
            Integer existingQual = qualArray[arrayPos][idx];
            if (existingQual < qual)
            {
                arr[arrayPos][idx] = base;
            }
            else if (existingQual == qual)
            {
                logger.warn("conflicting forward/reverse read bases: " + pi.getRecord().getReadName() + ", " + pi.getRefPosition() + ", " + arrayPos + ", " + idx + ", " + existing + ", " + base + ", " + qual);
                arr[arrayPos][idx] = 'X';
            }
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> ret = new ArrayList<>();
        ret.addAll(super.getDefaultReadFilters());
        ret.add(new ReadFilterLibrary.MappedReadFilter());

        if (minMappingQual > 0) {
            ret.add(new MappingQualityReadFilter(minMappingQual));
        }

        ret.add(new ReadFilterLibrary.PrimaryLineReadFilter());

        return Collections.unmodifiableList(ret);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Total read pairs inspected: " + readsInspected);
        logger.info("Total reads failing filters: " + failedFilters);
        logger.info("Total reads dropped for incomplete coverage: " + totalDroppedForCoverage);
        logger.info("Total reads by pairing:");
        readTotalHist.forEach((x, y) -> {
            logger.info(x + ": " + y);
        });

        for (SimpleInterval i : resultMap.keySet()) {
            outputStream.println("*******************************************");
            outputStream.println("Interval: " + i.toString());
            outputStream.println("Total read pairs inspected: " + readsInspected);
            outputStream.println("Total reads by reads/alignments per group:");
            readTotalHist.forEach((x, y) -> {
                outputStream.println('\t' + getPairLabel(x) + ": " + y);
            });

            outputStream.println("");
            try (IndexedFastaSequenceFile idx = new IndexedFastaSequenceFile(referenceArguments.getReferencePath()))
            {
                ReferenceSequence ref = idx.getSubsequenceAt(i.getContig(), i.getStart(), i.getEnd());
                Map<Character[][], Integer> haplotypes = resultMap.get(i);
                Map<Integer, TreeSet<Integer>> indels = getInsertionMap(haplotypes);
                String referenceSequence = getReferenceSequence(ref, indels);
                outputStream.println(referenceSequence);

                //convert to strings:
                Map<String, Integer> stringMap = new HashMap<>();
                for (Character[][] haplo : haplotypes.keySet()) {
                    String haplotypeSequence = convertHaplotypeToString(haplo, ref.getBases(), indels);
                    if (stringMap.containsKey(haplotypeSequence)){
                        throw new GATKException.ShouldNeverReachHereException("The map contains duplicate keys: " + haplotypeSequence);
                    }

                    stringMap.put(haplotypeSequence, haplotypes.get(haplo));
                }

                List<String> orderedList = new ArrayList<>(stringMap.keySet());
                Collections.sort(orderedList);

                orderedList.sort((a, b) -> {
                    return stringMap.get(b).compareTo(stringMap.get(a));
                });

                AtomicInteger totalHaplotypes = new AtomicInteger();
                haplotypes.forEach((x, y) -> {totalHaplotypes.addAndGet(y);});

                for (String haplotypeSequence : orderedList) {
                    outputStream.println(haplotypeSequence + '\t' + stringMap.get(haplotypeSequence) + '\t' + Utils.formattedPercent(stringMap.get(haplotypeSequence), totalHaplotypes.get()));
                }
            }
            catch (IOException e) {
                throw new GATKException("Unable to reference file: " + e.getMessage(), e);
            }

            outputStream.println("");
        }

        return super.onTraversalSuccess();
    }

    private String getPairLabel(Integer reads) {
        switch (reads) {
            case 1:
                return "Singleton";
            case 2:
                return "Paired";
            default:
                return reads.toString();
        }
    }

    @Override
    public void closeTool() {
        super.closeTool();

        if (bamReader != null){
            try {
                bamReader.close();
            }
            catch (IOException e) {
                //ignore
            }
        }

        outputStream.close();
    }

    private Map<Integer, TreeSet<Integer>> getInsertionMap(Map<Character[][], Integer> combinedResults)
    {
        //build list of all insertions that are present
        Map<Integer, TreeSet<Integer>> indels = new HashMap<>();
        for (Character[][] haplotype : combinedResults.keySet())
        {
            for (int idx = 0;idx < haplotype.length;idx++)
            {
                Character[] arr = haplotype[idx];
                if (arr == null)
                {
                    continue;
                }

                if (arr.length > 1)
                {
                    for (int i = 1; i < arr.length; i++)
                    {
                        TreeSet<Integer> l = indels.containsKey(idx) ? indels.get(idx) : new TreeSet<>();
                        l.add(i);
                        indels.put(idx, l);
                    }
                }
            }
        }

        return indels;
    }

    private String convertHaplotypeToString(Character[][] haplotype, byte[] refBases, Map<Integer, TreeSet<Integer>> indels)
    {
        StringBuilder sb = new StringBuilder();
        for (int idx = 0;idx < haplotype.length;idx++)
        {
            Character[] arr = haplotype[idx];
            char ref = Character.toUpperCase((char)refBases[idx]);

            if (arr == null)
            {
                sb.append(':');
            }
            else if (arr[0] == ref)
            {
                sb.append('.');
            }
            else
            {
                sb.append(arr[0]);
            }

            if (indels.containsKey(idx))
            {
                for (int insertIdx : indels.get(idx))
                {
                    if (arr != null && insertIdx < arr.length)
                    {
                        sb.append(arr[insertIdx]);
                    }
                    else
                    {
                        sb.append("-");
                    }
                }
            }
        }

        return sb.toString();
    }

    private String getReferenceSequence(ReferenceSequence ref, Map<Integer, TreeSet<Integer>> indels)
    {
        StringBuilder sb = new StringBuilder();
        int idx = 0;
        for (byte b : ref.getBases())
        {
            sb.append((char)b);
            if (indels.containsKey(idx))
            {
                for (int insertIdx : indels.get(idx))
                {
                    sb.append("-");
                }
            }

            idx++;
        }

        return sb.toString();
    }
}
