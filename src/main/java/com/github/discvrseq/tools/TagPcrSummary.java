package com.github.discvrseq.tools;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.util.SequenceMatcher;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.collections4.ComparatorUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.TextFeature;
import org.biojava.nbio.core.sequence.io.GenbankWriter;
import org.biojava.nbio.core.sequence.io.GenericGenbankHeaderFormat;
import org.biojava.nbio.core.sequence.location.SequenceLocation;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * This tool is designed to inspect the results of Tag-PCR (an assay to measure transgene integration into the genome); however, in theory it could be used with any assay
 * providing similar data.  It provides a very detailed output, but makes some assumptions and requires specific inputs:
 *
 * <ul>
 *     <li>The BAM is queryName sorted and the alignments per read inspected.  The tool assumes the genome is comprised of the organism plus one additional contig representing the delivery vector</li>
 *     <li>Only first-mate reads are considered.  Reverse are ignored.</li>
 *     <li>Each alignment is inspected (including cropped bases) for the presence of a short sequence expected to be at the insert/genome junction</li>
 *     <li>*IMPORTANT* The orientation of each hit is used to determine the orientation of the transgene in the genome.  Each query sequence is associated with one end of the transgene</li>
 *     <li>The hits are summarized based on the number of reads per genomic position (junction border)</li>
 *     <li>If --primer3-path and --primer-pair-table are provided, the tool will iterate each passing integration site, extract the upstream and downstream region (+/- 1000bp), and design primer pairs that site inside the transgene and flanking genomic region.</li>
 *     <li>If --blastn-path and --blast-db-path are provided, the tool will BLAST putative primers against the reference database and any primer with multiple hits will be flagged/discarded, along with non full-length primers that have a perfect match at the 3' end.</li>
 *     <li>To aid in inspecting the results, a genbank file can also be created (--genbank-output), which has one record per insert region, with the transgene region highlighted.  If primers were designed, these will also appear.</li>
 * </ul>
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar TagPcrSummary \
 *     -R currentGenome.fasta \
 *     -b myBam.bam \
 *     --output-table output.txt \
 *     --primer-pair-table primer_summary.txt \
 *     --prime3-path /usr/bin/primer3_core \
 *     --genbank-output output.gb
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "The is a specialist tool designed to summarize and interpret integration sites of a transgene into a genome",
        oneLineSummary = "Detect and summarize transgene integration",
        programGroup = DiscvrSeqDevProgramGroup.class
)
public class TagPcrSummary extends GATKTool {
    @Argument(fullName = "bam", shortName = "b", doc = "A BAM file with alignments to be inspected", common = false, optional = false)
    public File inputBam;

    @Argument(doc="File to which TSV output should be written", fullName = "output-table", shortName = "o", optional = false)
    public File outputTsv = null;

    @Argument(doc="The minimum number of alignments required at a position to report", fullName = "min-alignments", shortName = "ma", optional = true)
    public int MIN_ALIGNMENTS = 3;

    @Argument(doc="File to which a summary of integration sites, flanking sequence and optionally primers should be written in genbank format", fullName = "genbank-output", shortName = "g", optional = true)
    public File outputGenbank = null;

    @Argument(doc="File to which TSV summarizing potential primer pairs should be written", fullName = "primer-pair-table", shortName = "pt", optional = true)
    public File primerPairTable = null;

    @Argument(doc="In order for this tool to design validation primers, the path to primer3 must be provided", fullName = "primer3-path", shortName = "p3", optional = true)
    public String primer3ExePath = "primer3";

    @Argument(doc="File to which a TSV of summary metrics should be written", fullName = "metrics-table", shortName = "mt", optional = true)
    public File metricsFile = null;

    @Argument(doc="The path to the blastn executable.  This is required for BLAST validation to be perform against putative primers.", fullName = "blastn-path", shortName = "bn", optional = true)
    public String blastnPath = "blastn";

    @Argument(doc="In order for this tool to use BLAST to detect validate the primers by detecting alternate binding sites, the path to a BLAST DB compiled against this reference FASTA must be provided", fullName = "blast-db-path", shortName = "bdb", optional = true)
    public String blastDatabase = null;

    @Argument(doc="If BLAST will be used, this value is passed to the -num_threads argument of blastn.", fullName = "blast-threads", shortName = "bt", optional = true)
    public Integer blastThreads = null;

    @Override
    public boolean requiresReference() {
        return true;
    }

    List<InsertDescriptor> INSERT_DESCRIPTORS = new ArrayList<>();

    @Override
    protected void onStartup() {
        super.onStartup();

        INSERT_DESCRIPTORS.add(new InsertDescriptor("pENTR-PB511-Repro", "pENTR-PB511-Repro", Arrays.asList(
                new InsertJunctionDescriptor("PB-3TR", Collections.singletonList("GCAGACTATCTTTCTAGGGTTAA"), new Interval("pENTR-PB511-Repro", 7345, 7555), "PB-3TR", new Interval("pENTR-PB511-Repro", 606, 1066), "PB-5TR", Collections.singletonList("CATTGACAAGCACGCCTCAC")),
                new InsertJunctionDescriptor("PB-5TR", Collections.singletonList("ATGATTATCTTTCTAGGGTTAA"), new Interval("pENTR-PB511-Repro", 606, 1066), "PB-5TR", new Interval("pENTR-PB511-Repro", 7345, 7555), "PB-3TR", Arrays.asList("CACATGATTATCTTTAACGTACGTCAC", "GACCGATAAAACACATGCGTCA"))
        )));

        if (blastDatabase != null) {
            File blastTest = new File(blastDatabase + ".nhr");
            if (!blastTest.exists()) {
                throw new UserException.BadInput("Invalid BLAST DB path.  Expected this location to contain many files, including: " + blastTest.getPath());
            }
        }
    }

    @Override
    public void traverse() {
        NumberFormat pf = NumberFormat.getPercentInstance();
        pf.setMaximumFractionDigits(6);

        SamReaderFactory fact = SamReaderFactory.makeDefault();
        fact.validationStringency(ValidationStringency.SILENT);
        fact.referenceSequence(referenceArguments.getReferencePath());

        int totalAlignments = 0;
        int totalReads = 0;
        int numReadsSpanningJunction = 0;
        File bam = ensureQuerySorted(inputBam);

        Map<String, JunctionMatch> totalMatches = new HashMap<>();
        Map<String, Number> metricsMap = new HashMap<>();

        Map<String, Integer> primaryAlignmentsToInsert = new HashMap<>();
        Map<String, Integer> secondaryAlignmentsToInsert = new HashMap<>();
        int reverseReadsSkipped = 0;

        try (SamReader bamReader = fact.open(bam); SAMRecordIterator it = bamReader.iterator()) {

            List<SAMRecord> alignmentsForRead = new ArrayList<>();
            while (it.hasNext()) {
                SAMRecord rec = it.next();
                if (rec.getReadUnmappedFlag()) {
                    continue;
                }

                // Skip reverse reads
                if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                    reverseReadsSkipped++;
                    continue;
                }

                totalAlignments++;
                if (alignmentsForRead.isEmpty() || alignmentsForRead.get(0).getReadName().equals(rec.getReadName())) {
                    alignmentsForRead.add(rec);
                } else {
                    totalReads++;
                    numReadsSpanningJunction += processAlignmentsForRead(alignmentsForRead, totalMatches, primaryAlignmentsToInsert, secondaryAlignmentsToInsert);
                }
            }

            //ensure we capture final read
            if (!alignmentsForRead.isEmpty()) {
                totalReads++;
                numReadsSpanningJunction += processAlignmentsForRead(alignmentsForRead, totalMatches, primaryAlignmentsToInsert, secondaryAlignmentsToInsert);

            }
        }
        catch (IOException e)
        {
            throw new GATKException(e.getMessage(), e);
        }

        double pct = totalReads == 0 ? 0 : numReadsSpanningJunction / (double)totalReads;
        logger.info("Total reads spanning a junction: " + numReadsSpanningJunction + " of " + totalReads + " (" + pf.format(pct) + ")");

        metricsMap.put("TotalAlignments", totalAlignments);
        metricsMap.put("DistinctReads", totalReads);
        metricsMap.put("NumReadsSpanningJunction", numReadsSpanningJunction);
        metricsMap.put("PctReadsSpanningJunction", pct);

        Set<String> inserts = new TreeSet<>();
        inserts.addAll(primaryAlignmentsToInsert.keySet());
        inserts.addAll(secondaryAlignmentsToInsert.keySet());
        logger.info("Total primary/secondary alignments matching inserts:");
        int totalPrimaryAlignmentsMatchingInsert = 0;
        int totalSecondaryAlignmentsMatchingInsert = 0;
        for (String key : inserts) {
            logger.info(key + ", primary: " + primaryAlignmentsToInsert.getOrDefault(key, 0));
            totalPrimaryAlignmentsMatchingInsert += primaryAlignmentsToInsert.getOrDefault(key, 0);

            logger.info(key + ", secondary: " + secondaryAlignmentsToInsert.getOrDefault(key, 0));
            totalPrimaryAlignmentsMatchingInsert += secondaryAlignmentsToInsert.getOrDefault(key, 0);
        }

        metricsMap.put("TotalPrimaryAlignmentsMatchingInsert", totalPrimaryAlignmentsMatchingInsert);
        metricsMap.put("FractionPrimaryAlignmentsMatchingInsert", (double)totalPrimaryAlignmentsMatchingInsert / totalReads);
        metricsMap.put("TotalSecondaryAlignmentsMatchingInsert", totalSecondaryAlignmentsMatchingInsert);

        int totalReverse = 0;
        List<DNASequence> amplicons = new ArrayList<>();

        if (!totalMatches.isEmpty()) {
            try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(outputTsv), '\t', CSVWriter.NO_QUOTE_CHARACTER); CSVWriter primerWriter = primerPairTable == null ? null : new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(primerPairTable), '\t', CSVWriter.NO_QUOTE_CHARACTER); ReferenceSequenceFile refSeq = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceArguments.getReferencePath())) {
                writer.writeNext(new String[]{"SiteName", "JunctionName", "Orientation", "Chr", "Position", "Orientation", "Total"});
                if (primerWriter != null) {
                    primerWriter.writeNext(new String[]{"SiteName", "JunctionName", "Orientation", "Primer-F-Name", "Primer-F-Seq", "Primer-F-Start", "Primer-R-Name", "Primer-R-Seq", "Primer-R-Start"});
                }
                Map<String, ReferenceSequence> refMap = new HashMap<>();
                List<String> sortedKeys = new ArrayList<>(totalMatches.keySet());
                sortedKeys.sort(ComparatorUtils.naturalComparator());

                int totalOutput = 0;
                for (String key : sortedKeys) {
                    JunctionMatch jm = totalMatches.get(key);

                    if (jm.totalReads >= MIN_ALIGNMENTS) {
                        totalOutput++;

                        ReferenceSequence ref = refMap.get(jm.contigName);
                        if (ref == null) {
                            ref = refSeq.getSequence(jm.contigName);
                            refMap.put(jm.contigName, ref);
                        }

                        if (jm.matchIsNegativeStrand) {
                            totalReverse++;
                        }

                        String siteName = "Site" + totalOutput;
                        writer.writeNext(new String[]{siteName, jm.jd.junctionName, (jm.matchIsNegativeStrand ? "Minus" : "Plus"), jm.contigName, String.valueOf(jm.alignmentStart), (jm.matchIsNegativeStrand ? "-" : "+"), String.valueOf(jm.totalReads)});

                        Pair<DNASequence, DNASequence> ampliconPair = jm.getAmplicons(refSeq, refMap, siteName, primerWriter);
                        amplicons.add(ampliconPair.getKey());
                        amplicons.add(ampliconPair.getValue());

                    }
                }

                logger.info("Total junctions output: " + totalOutput);
                logger.info("Forward / Reverse: " + (totalOutput - totalReverse) + " / " + totalReverse);

                metricsMap.put("TotalIntegrationSitesOutput", totalOutput);
                metricsMap.put("IntegrationSitesOutputMinusStrand", totalReverse);
                metricsMap.put("IntegrationSitesOutputPlusStrand", (totalOutput - totalReverse));
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

        }

        if (blastDatabase != null && primerPairTable != null) {
            runBlastN(primerPairTable, amplicons);
        }

            //Now make genbank output:
            if (outputGenbank != null) {
                try (BufferedOutputStream os = new BufferedOutputStream(IOUtil.openFileForWriting(outputGenbank))) {
                    GenbankWriter<DNASequence, NucleotideCompound> genbankWriter = new GenbankWriter<>(os, amplicons, new GenericGenbankHeaderFormat<>("LINEAR"));
                    genbankWriter.process();
                }
                catch (Exception e) {
                    throw new GATKException(e.getMessage(), e);
                }
            }


        if (metricsFile != null && !metricsMap.isEmpty()) {
            try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(metricsFile), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
                writer.writeNext(new String[]{"MetricName", "MetricValue"});
                metricsMap.forEach((key, value) -> writer.writeNext(new String[]{key, String.valueOf(value)}));
            }
            catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }
        }
    }

    private int processAlignmentsForRead(List<SAMRecord> alignmentsForRead, Map<String, JunctionMatch> totalMatches, Map<String, Integer> primaryAlignmentsToInsert, Map<String, Integer> secondaryAlignmentsToInsert) {
        Map<String, JunctionMatch> matches = new HashMap<>();
        alignmentsForRead.forEach(rec -> {
            if (rec.isSecondaryOrSupplementary()) {
                INSERT_DESCRIPTORS.forEach(id -> {
                    if (rec.getContig().equals(id.contigName)) {
                        int total = secondaryAlignmentsToInsert.getOrDefault(id.displayName, 0);
                        total++;
                        secondaryAlignmentsToInsert.put(id.displayName, total);
                    }
                });
            }
            else {
                for (InsertDescriptor id : INSERT_DESCRIPTORS) {
                    for (InsertJunctionDescriptor jd : id.junctions) {
                        matches.putAll(jd.getMatches(rec, id));
                    }
                }

                INSERT_DESCRIPTORS.forEach(id -> {
                    if (rec.getContig().equals(id.contigName)) {
                        int total = primaryAlignmentsToInsert.getOrDefault(id.displayName, 0);
                        total++;
                        primaryAlignmentsToInsert.put(id.displayName, total);
                    }
                });
            }
        });

        alignmentsForRead.clear();
        if (!matches.isEmpty()) {
            matches.forEach((x, y) -> {
                if (totalMatches.containsKey(x)) {
                    totalMatches.get(x).addRead();
                }
                else {
                    totalMatches.put(x, y);
                }
            });

            return 1;
        } else {
            return 0;
        }
    }

    private File ensureQuerySorted(File bam) {
        SamReaderFactory fact = SamReaderFactory.makeDefault();
        fact.validationStringency(ValidationStringency.SILENT);
        fact.referenceSequence(referenceArguments.getReferencePath());

        SAMFileHeader.SortOrder so = fact.getFileHeader(bam).getSortOrder();
        if (so == SAMFileHeader.SortOrder.queryname) {
            logger.info("BAM is already query sorted, no need to sort");
            return bam;
        }

        logger.info("Sorting BAM in queryName order");
        try (SamReader reader = SamReaderFactory.makeDefault().referenceSequence(referenceArguments.getReferencePath()).open(bam)) {
            reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.queryname);

            File querySorted = File.createTempFile(bam.getName(), ".querySorted.bam").toPath().normalize().toFile();

            try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, querySorted)) {
                for (final SAMRecord rec : reader) {
                    writer.addAlignment(rec);
                }
            }

            return querySorted;
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }
    }

    public class InsertDescriptor {
        final String displayName;
        final String contigName;

        final List<InsertJunctionDescriptor> junctions;

        public InsertDescriptor(String displayName, String contigName, List<InsertJunctionDescriptor> junctions) {
            this.displayName = displayName;
            this.contigName = contigName;
            this.junctions = junctions;
        }
    }

    public class InsertJunctionDescriptor {
        final String junctionName;
        final List<String> searchStrings;
        List<String> searchStringsRC = null;
        int mismatchesAllowed = 2;

        final Interval insertRegionProximal;
        final String insertRegionProximalLabel;

        final Interval insertRegionDistal;
        String insertRegionDistalLabel;
        final List<String> internalPrimers;

        public InsertJunctionDescriptor(String junctionName, List<String> searchStrings, Interval insertRegionProximal, String insertRegionProximalLabel, Interval insertRegionDistal, String insertRegionDistalLabel, List<String> internalPrimers) {
            this.junctionName = junctionName;
            this.searchStrings = Collections.unmodifiableList(searchStrings);
            this.insertRegionProximal = insertRegionProximal;
            this.insertRegionProximalLabel = insertRegionProximalLabel;
            this.insertRegionDistal = insertRegionDistal;
            this.insertRegionDistalLabel = insertRegionDistalLabel;
            this.internalPrimers = Collections.unmodifiableList(internalPrimers);
        }

        public String getFeatureLabel(END_TYPE type) {
            switch (type) {
                case promimal:
                    return insertRegionProximalLabel;
                case distal:
                    return insertRegionDistalLabel;
                default:
                    throw new IllegalArgumentException("Unknown type");
            }
        }

        public void setMismatchesAllowed(int mismatchesAllowed) {
            this.mismatchesAllowed = mismatchesAllowed;
        }

        private List<String> getSearchStrings(boolean reverseComplement) {
            if (reverseComplement) {
                if (searchStringsRC == null) {
                    searchStringsRC = new ArrayList<>();
                    searchStrings.forEach(x -> {
                        searchStringsRC.add(SequenceUtil.reverseComplement(x));
                    });
                }

                return searchStringsRC;
            }

            return searchStrings;
        }

        public Map<String, JunctionMatch> getMatches(SAMRecord rec, InsertDescriptor id) {
            Map<String, JunctionMatch> matches = new HashMap<>();

            //Forward orientation.
            for (String query : getSearchStrings(false)) {
                Integer match = SequenceMatcher.fuzzyMatch(query, rec.getReadString().toUpperCase(), mismatchesAllowed);
                if (match != null) {
                    //let the read deal with clipping
                    JunctionMatch jm = new JunctionMatch(this, rec.getContig(), rec.getAlignmentStart(), false, rec.getContig().equals(id.contigName));
                    matches.put(jm.getKey(), jm);
                }
            }

            //Reverse:
            for (String query : getSearchStrings(rec.getReadNegativeStrandFlag())) {
                Integer match = SequenceMatcher.fuzzyMatch(query, rec.getReadString().toUpperCase(), mismatchesAllowed);
                if (match != null) {
                    JunctionMatch jm = new JunctionMatch(this, rec.getContig(), rec.getAlignmentEnd(), true, rec.getContig().equals(id.contigName));
                    matches.put(jm.getKey(), jm);
                }
            }

            return matches;
        }
    }

    private enum END_TYPE {
        distal(),
        promimal()
    }

    private class JunctionMatch {
        private InsertJunctionDescriptor jd;

        private String contigName;
        private int alignmentStart;
        private boolean matchIsNegativeStrand;
        private boolean matchesInsertContig;
        private int totalReads = 1;

        public JunctionMatch(InsertJunctionDescriptor jd, String contigName, int alignmentStart, boolean matchIsNegativeStrand, boolean matchesInsertContig) {
            this.jd = jd;
            this.contigName = contigName;
            this.alignmentStart = alignmentStart;
            this.matchIsNegativeStrand = matchIsNegativeStrand;
            this.matchesInsertContig = matchesInsertContig;
        }

        public String getKey() {
            return jd.junctionName + "<>" + contigName + "<>" + alignmentStart + "<>" + matchIsNegativeStrand;
        }

        public void addRead() {
            totalReads += 1;
        }

        public Pair<DNASequence, DNASequence> getAmplicons(ReferenceSequenceFile refSeq, Map<String, ReferenceSequence> refMap, String siteName, CSVWriter primerWriter) {
            ReferenceSequence ref = refMap.get(contigName);
            if (ref == null) {
                ref = refSeq.getSequence(contigName);
                refMap.put(contigName, ref);
            }

            Interval afterInterval = new Interval(ref.getName(), alignmentStart, Math.min(alignmentStart + 1000, ref.length()));
            String after = ref.getBaseString().substring(afterInterval.getStart() - 1, afterInterval.getEnd());

            Interval beforeInterval = new Interval(ref.getName(), Math.max(1, alignmentStart - 1000), alignmentStart);
            String before = ref.getBaseString().substring(beforeInterval.getStart() - 1, beforeInterval.getEnd());

            String insertRegionProximal = "";
            if (jd.insertRegionProximal != null) {
                ReferenceSequence insertRef = refMap.get(jd.insertRegionProximal.getContig());
                if (insertRef == null) {
                    insertRef = refSeq.getSequence(jd.insertRegionProximal.getContig());
                    refMap.put(insertRef.getName(), insertRef);
                }

                insertRegionProximal = insertRef.getBaseString().substring(jd.insertRegionProximal.getStart() - 1, jd.insertRegionProximal.getEnd());
            }

            String insertRegionDistal = "";
            if (jd.insertRegionDistal != null) {
                ReferenceSequence insertRef = refMap.get(jd.insertRegionDistal.getContig());
                if (insertRef == null) {
                    insertRef = refSeq.getSequence(jd.insertRegionDistal.getContig());
                    refMap.put(insertRef.getName(), insertRef);
                }

                insertRegionDistal = insertRef.getBaseString().substring(jd.insertRegionDistal.getStart() - 1, jd.insertRegionDistal.getEnd());
            }

            DNASequence seq1;
            DNASequence seq2;
            try {
                String id = siteName;
                String idFull = contigName + ":" + alignmentStart + ":" + (matchIsNegativeStrand ? "Minus" : "Plus");

                //Note: this is a limitation in biojava
                if (id.length() > 16) {
                    logger.error("A site name greater than 16 characters was passed: " + id + ".  This will be truncated to comply with genbank format");
                    id = id.substring(0, 16);
                }

                String orientationSuffix = matchIsNegativeStrand ? "m" : "";
                if (matchIsNegativeStrand) {
                    seq1 = new DNASequence(SequenceUtil.reverseComplement(insertRegionDistal) + after);
                    String ampliconName = id + "-" + jd.getFeatureLabel(END_TYPE.distal);
                    String accession = ampliconName + orientationSuffix;
                    if (accession.length() > 15) {
                        //a limitation of BioJava
                        accession = accession.substring(0, 15);
                    }
                    seq1.setAccession(new AccessionID(accession));
                    seq1.setDescription(idFull);

                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(jd.getFeatureLabel(END_TYPE.distal) + "-RC", "Vector", jd.junctionName, jd.insertRegionDistal.toString());
                    tf.setLocation(new SequenceLocation<>(1, insertRegionDistal.length(), seq1, Strand.NEGATIVE));
                    seq1.addFeature(tf);
                    runPrimer3(siteName, jd.getFeatureLabel(END_TYPE.distal), seq1, insertRegionDistal.length(), outputTsv.getParentFile(), primerWriter);

                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tfG = new TextFeature<>(afterInterval.getContig() + ":" + afterInterval.getStart(), "Genome", afterInterval.toString(), afterInterval.toString());
                    tfG.setLocation(new SequenceLocation<>(insertRegionDistal.length() + 1, seq1.getBioEnd(), seq1, Strand.POSITIVE));
                    seq1.addFeature(tfG);

                    seq2 = new DNASequence(before + SequenceUtil.reverseComplement(insertRegionProximal));
                    String ampliconName2 = id + "-" + jd.getFeatureLabel(END_TYPE.promimal);
                    String accession2 = ampliconName2 + orientationSuffix;
                    if (accession2.length() > 15) {
                        //a limitation of BioJava
                        accession2 = accession2.substring(0, 15);
                    }
                    seq2.setAccession(new AccessionID(accession2));
                    seq2.setDescription(idFull);

                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf2 = new TextFeature<>(jd.getFeatureLabel(END_TYPE.promimal) + "-RC", "Vector", jd.junctionName, jd.insertRegionProximal.toString());
                    tf2.setLocation(new SequenceLocation<>(before.length() + 1, seq2.getBioEnd(), seq2, Strand.NEGATIVE));
                    seq2.addFeature(tf2);
                    runPrimer3(siteName, jd.getFeatureLabel(END_TYPE.promimal), seq2, before.length(), outputTsv.getParentFile(), primerWriter);

                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tfG2 = new TextFeature<>(beforeInterval.getContig() + ":" + beforeInterval.getStart(), "Genome", beforeInterval.toString(), beforeInterval.toString());
                    tfG2.setLocation(new SequenceLocation<>(1, before.length(), seq2, Strand.POSITIVE));
                    seq2.addFeature(tfG2);
                } else {
                    seq1 = new DNASequence(before + insertRegionProximal);
                    String ampliconName = id + "-" + jd.getFeatureLabel(END_TYPE.promimal);
                    String accession = ampliconName + orientationSuffix;
                    if (accession.length() > 15) {
                        //a limitation of BioJava
                        accession = accession.substring(0, 15);
                    }
                    seq1.setAccession(new AccessionID(accession));
                    seq1.setDescription(idFull);
                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(jd.getFeatureLabel(END_TYPE.promimal), "Vector", jd.junctionName, jd.insertRegionProximal.toString());
                    tf.setLocation(new SequenceLocation<>(before.length() + 1, seq1.getBioEnd(), seq1, Strand.POSITIVE));
                    seq1.addFeature(tf);
                    runPrimer3(siteName, jd.getFeatureLabel(END_TYPE.promimal), seq1, before.length(), outputTsv.getParentFile(), primerWriter);

                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tfG = new TextFeature<>(afterInterval.getContig() + ":" + afterInterval.getStart(), "Genome", afterInterval.toString(), afterInterval.toString());
                    tfG.setLocation(new SequenceLocation<>(1, before.length(), seq1, Strand.POSITIVE));
                    seq1.addFeature(tfG);

                    seq2 = new DNASequence(insertRegionDistal + after);
                    String ampliconName2 = id + "-" + jd.getFeatureLabel(END_TYPE.distal);
                    String accession2 = ampliconName2 + orientationSuffix;
                    if (accession2.length() > 15) {
                        //a limitation of BioJava
                        accession2 = accession2.substring(0, 15);
                    }
                    seq2.setAccession(new AccessionID(accession2));
                    seq2.setDescription(idFull);
                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf2 = new TextFeature<>(jd.getFeatureLabel(END_TYPE.distal), "Vector", jd.junctionName, jd.insertRegionDistal.toString());
                    tf2.setLocation(new SequenceLocation<>(1, insertRegionDistal.length(), seq2, Strand.POSITIVE));
                    seq2.addFeature(tf2);
                    runPrimer3(siteName, jd.getFeatureLabel(END_TYPE.distal), seq2, insertRegionDistal.length(), outputTsv.getParentFile(), primerWriter);

                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tfG2 = new TextFeature<>(beforeInterval.getContig() + ":" + beforeInterval.getStart(), "Genome", beforeInterval.toString(), beforeInterval.toString());
                    tfG2.setLocation(new SequenceLocation<>(insertRegionDistal.length() + 1, seq2.getBioEnd(), seq2, Strand.POSITIVE));
                    seq2.addFeature(tfG2);
                }
            } catch (CompoundNotFoundException e) {
                throw new GATKException(e.getMessage(), e);
            }

            return Pair.of(seq1, seq2);
        }

        private void runPrimer3(String siteName, String junctionName, DNASequence sequence, int junctionSite, File outDir, CSVWriter primerWriter) {
            logger.info("Running primer3 for: " + siteName + ", " + junctionName);

            String seqString = sequence.getSequenceAsString();

            List<String> args = new ArrayList<>();
            args.add(primer3ExePath);

            File output = new File(outDir, siteName + "-" + junctionName + ".p3.output.txt");
            args.add("--output=" + output.getPath());

            args.add("--default_version=2");
            args.add("--strict_tags");

            File errorFile = new File(outDir, "p3.error.txt");
            args.add("--error=" + errorFile.getPath());

            int start = junctionSite - 50;
            File input = new File(outDir, siteName + "-" + junctionName + ".p3.input.txt");
            try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(input)) {
                writer.write(
                        //TODO: SEQUENCE_EXCLUDED_REGION 10,20 50,10
                        "SEQUENCE_ID=" + sequence.getAccession().getID() + "\n" +
                                "SEQUENCE_TEMPLATE=" + seqString + "\n" +
                                "SEQUENCE_TARGET=" + start + ",100" + "\n" +
                                "PRIMER_PICK_LEFT_PRIMER=1" + "\n" +
                                "PRIMER_PICK_INTERNAL_OLIGO=0" + "\n" +
                                "PRIMER_PICK_RIGHT_PRIMER=1" + "\n" +
                                "PRIMER_OPT_SIZE=18" + "\n" +
                                "PRIMER_MIN_SIZE=15" + "\n" +
                                "PRIMER_MAX_SIZE=21" + "\n" +
                                "PRIMER_MAX_NS_ACCEPTED=1" + "\n" +
                                "PRIMER_PRODUCT_SIZE_RANGE=100-300 501-600 601-700 401-500 701-850" + "\n" +
                                "P3_FILE_FLAG=0" + "\n" +
                                "PRIMER_EXPLAIN_FLAG=1" + "\n" +
                                "="
                );
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            args.add(input.getPath());

            final ProcessSettings prs = new ProcessSettings(args.toArray(new String[args.size()]));
            prs.getStderrSettings().printStandard(true);
            prs.getStdoutSettings().printStandard(true);
            ProcessOutput po = ProcessController.getThreadLocal().exec(prs);
            if (po.getExitValue() != 0) {
                logger.info("primer3 had a non-zero exit");
            }

            if (!output.exists()) {
                throw new GATKException("Error running primer3");
            }

            Map<String, String> outputMap = new HashMap<>();
            try (BufferedReader reader = IOUtil.openFileForBufferedUtf8Reading(output)) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.isEmpty() || "=".equals(line)) {
                        continue;
                    }

                    String[] tokens = line.split("=");
                    outputMap.put(tokens[0], tokens[1]);
                }
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            input.delete();
            errorFile.delete();

            if (!outputMap.containsKey("PRIMER_ERROR")) {
                logger.info("primer3 error: " + outputMap.get("PRIMER_ERROR"));
            }

            if (!outputMap.containsKey("PRIMER_PAIR_NUM_RETURNED")) {
                logger.info("primer3 returned no pairs");
                return;
            }

            output.delete();

            int numPairs = Integer.parseInt(outputMap.get("PRIMER_PAIR_NUM_RETURNED"));
            logger.info("total candidate pairs: " + numPairs);

            for (int i = 0; i < numPairs; i++) {
                String primer1Seq = outputMap.get("PRIMER_LEFT_" + i + "_SEQUENCE");
                int primer1Loc = Integer.parseInt(outputMap.get("PRIMER_LEFT_" + i).split(",")[0]) + 1;

                String primer2Seq = outputMap.get("PRIMER_RIGHT_" + i + "_SEQUENCE");
                int primer2Loc = Integer.parseInt(outputMap.get("PRIMER_RIGHT_" + i).split(",")[0]) + 1;

                String pairNameShort = junctionName + "-" + (i + 1);
                String pairName = siteName + "-" + pairNameShort;
                if (primerWriter != null) {
                    primerWriter.writeNext(new String[]{siteName, junctionName, (matchIsNegativeStrand ? "Minus" : "Plus"), pairName + "F", primer1Seq, String.valueOf(primer1Loc), pairName + "R", primer2Seq, String.valueOf(primer2Loc)});
                }

                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(pairNameShort + "F", "Primer3", pairName + "F", pairName + "F");
                tf.setLocation(new SequenceLocation<>(primer1Loc, primer1Loc + primer1Seq.length(), sequence, Strand.POSITIVE));
                sequence.addFeature(tf);

                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf2 = new TextFeature<>(pairNameShort + "R", "Primer3", pairName + "R", pairName + "R");
                tf2.setLocation(new SequenceLocation<>(primer2Loc, primer2Loc + primer2Seq.length(), sequence, Strand.NEGATIVE));
                sequence.addFeature(tf2);
            }
        }
    }

    private void runBlastN(File primer3Table, List<DNASequence> amplicons) {
        logger.info("Running BLAST to validate primers");
        File outDir = primer3Table.getParentFile();

        File blastInput = new File(outDir, "blastInput.fasta");
        File blastOutput = new File(outDir, "blastOutput.txt");

        List<String[]> primerLines = new ArrayList<>();
        Map<String, String> idxToPrimer = new HashMap<>();
        int primersSkipped = 0;
        int idx = 0;
        try (CSVReader reader = new CSVReader(IOUtil.openFileForBufferedUtf8Reading(primer3Table), '\t');PrintWriter writer = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(blastInput))) {
            String[] line;
            Set<String> primerSeqs = new HashSet<>();

            while ((line = reader.readNext()) != null) {
                if ("SiteName".equals(line[0])){
                    continue;
                }

                idx++;
                primerLines.add(line);

                if (!primerSeqs.contains(line[4].toUpperCase())) {
                    String name = "PrimerF-" + idx;
                    writer.println(">" + name);
                    writer.println(line[4].toUpperCase());
                    primerSeqs.add(line[4].toUpperCase());

                    idxToPrimer.put(name, line[4].toUpperCase());
                }
                else {
                    primersSkipped++;
                }

                if (!primerSeqs.contains(line[7].toUpperCase())) {
                    String name = "PrimerR-" + idx;
                    writer.println(">" + name);
                    writer.println(line[7].toUpperCase());
                    primerSeqs.add(line[7].toUpperCase());

                    idxToPrimer.put(name, line[7].toUpperCase());
                }
                else {
                    primersSkipped++;
                }
            }

            logger.info("Duplicate primers collapsed before BLAST: " + primersSkipped + " of " + (idx *2));
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        if (idxToPrimer.isEmpty()) {
            logger.info("No primers returned, skipping BLAST");
            return;
        }

        List<String> args = new ArrayList<>();
        args.add(blastnPath);

        args.add("-task");
        args.add("blastn-short");

        args.add("-db");
        args.add(blastDatabase);

        args.add("-query");
        args.add(blastInput.getPath());

        args.add("-out");
        args.add(blastOutput.getPath());

        args.add("-evalue");
        args.add("1");

        args.add("-max_hsps");
        args.add("2");

        args.add("-num_alignments");
        args.add("5");

        if (blastThreads != null) {
            args.add("-num_threads");
            args.add(blastThreads.toString());
        }

        args.add("-outfmt");
        args.add("6 qseqid sallseqid qstart qend sstart send evalue bitscore length nident sstrand");

        final ProcessSettings prs = new ProcessSettings(args.toArray(new String[args.size()]));
        prs.getStderrSettings().printStandard(true);
        prs.getStdoutSettings().printStandard(true);
        ProcessController.getThreadLocal().exec(prs);

        if (!blastOutput.exists()) {
            throw new GATKException("BLASTn did not produce an output, expected: " + blastOutput.getPath());
        }

        File blastOutputTable = new File(outDir, "blastResults.txt");
        Map<String, Boolean> failedPrimers = new HashMap<>();
        try (CSVReader reader = new CSVReader(IOUtil.openFileForBufferedUtf8Reading(blastOutput), '\t');CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(blastOutputTable), '\t', CSVWriter.NO_QUOTE_CHARACTER);) {
            String[] line;

            writer.writeNext(new String[]{"PrimerName", "Hit", "Start", "End", "Strand", "QueryStart", "QueryEnd", "EValue", "BitScore", "HitLength", "NumIdents", "FullLengthHit", "EndMatch"});
            Map<String, List<String[]>> hitMap = new LinkedHashMap<>();
            while ((line = reader.readNext()) != null) {
                String query = line[0];
                List<String[]> hits = hitMap.getOrDefault(query, new ArrayList<>());
                hits.add(line);
                hitMap.put(query, hits);
            }

            //Iterate BLAST hits, summarize by pair:
            hitMap.forEach((query, hit) -> {
                final String primerSeq = idxToPrimer.get(query);

                AtomicInteger endPrimeMatches = new AtomicInteger(0);
                List<String[]> perfectHits = new ArrayList<>();
                hit.forEach(h -> {
                    boolean fullLength = primerSeq.length() == Integer.parseInt(h[8]);
                    boolean endPrimeMatch = Integer.parseInt(h[3]) == primerSeq.length();
                    if (fullLength) {
                        perfectHits.add(h);
                    }
                    else if (endPrimeMatch) {
                        endPrimeMatches.getAndIncrement();
                    }

                    writer.writeNext(new String[]{primerSeq, h[1], h[4], h[5], h[10], h[2], h[3], h[6], h[7], h[8], h[9], String.valueOf(fullLength), String.valueOf(endPrimeMatch)});
                });

                if (perfectHits.size() > 1) {
                    failedPrimers.put(primerSeq, false);
                }
                else if (endPrimeMatches.get() > 0) {
                    failedPrimers.put(primerSeq, false);
                }
            });

            blastOutput.delete();
            blastInput.delete();
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        //Now write updated output:
        Set<String> primersToRemove = new HashSet<>();
        try (CSVWriter primerWriter = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(primer3Table), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            primerWriter.writeNext(new String[]{"SiteName", "JunctionName", "Orientation", "Primer-F-Name", "Primer-F-Seq", "Primer-F-Start", "Primer-R-Name", "Primer-R-Seq", "Primer-R-Start", "PassedBlast"});
            for (String[] line : primerLines) {
                boolean forwardPass = !failedPrimers.containsKey(line[4].toUpperCase());
                boolean reversePass = !failedPrimers.containsKey(line[7].toUpperCase());
                List<String> toWrite = new ArrayList<>(Arrays.asList(line));
                toWrite.add((forwardPass && reversePass) ? "Y" : "");
                primerWriter.writeNext(toWrite.toArray(new String[toWrite.size()]));

                if (!forwardPass || !reversePass) {
                    primersToRemove.add(line[3]);
                    primersToRemove.add(line[6]);
                }
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        logger.info("Total primer paired discarded due to multiple BLAST hits: " + (primersToRemove.size() / 2) + " of " + primerLines.size());

        //prune failed from genbank output:
        if (!primersToRemove.isEmpty()) {
            amplicons.listIterator().forEachRemaining(seq -> {
                List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = new ArrayList<>(seq.getFeatures());
                features.listIterator().forEachRemaining(feature -> {
                    if (primersToRemove.contains(feature.getDescription())) {
                        feature.setType("Fail:" + feature.getType());
                    }
                });
            });
        }
    }
}


