package com.github.discvrseq.walkers.tagpcr;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.github.discvrseq.tools.DiscvrSeqDevProgramGroup;
import com.github.discvrseq.util.SequenceMatcher;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.collections4.ComparatorUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.SystemUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.Logger;
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
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * This tool is designed to inspect the results of Tag-PCR (an assay to measure transgene integration into the genome); however, in theory it could be used with any assay
 * providing similar data.  It provides a very detailed output, but makes some assumptions and requires specific inputs:
 *
 * <ul>
 *     <li>The BAM is queryName sorted and the alignments per read inspected.  The tool assumes the genome is the reference for that organism. It is possible for this to also contain the transgene and/or delivery vector but this is not necessarily needed and can be more complicated to maintain</li>
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
 *
 * To get the complete output, the tool requires a file with a detailed description of the transgene and expected junction sites. Depending on your delivery system, this is likely specific to your plasmid/vector. Below are two examples:
 *
 * <pre>
 *
 *
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "The is a specialist tool designed to summarize and interpret integration sites of a transgene into a genome",
        oneLineSummary = "Detect and summarize transgene integration",
        programGroup = DiscvrSeqDevProgramGroup.class
)
public class TagPcrSummary extends GATKTool {
    @Argument(fullName = "bam", shortName = "b", doc = "A BAM file with alignments to be inspected", common = false, optional = true)
    public File inputBam;

    @Argument(doc="File to which TSV output should be written", fullName = "output-table", shortName = "o", optional = true)
    public File outputTsv = null;

    @Argument(doc="The minimum number of alignments required at a position to report", fullName = "min-alignments", shortName = "ma", optional = true)
    public int MIN_ALIGNMENTS = 3;

    @Argument(doc="Only sites with at least this fraction of total reads will be reported", fullName = "min-fraction", shortName = "mf", optional = true)
    public double MIN_FRACTION = 0.0;

    @Argument(doc="If greater than zero, up to this many reads will be written as a FASTA file for each site.  This can be useful to validate the junction border", fullName = "reads-to-output", shortName = "ro", optional = true)
    public int READS_PER_SITE = 0;

    @Argument(doc="File to which a summary of integration sites, flanking sequence and optionally primers should be written in genbank format", fullName = "genbank-output", shortName = "g", optional = true)
    public File outputGenbank = null;

    @Argument(doc="File to which TSV summarizing potential primer pairs should be written", fullName = "primer-pair-table", shortName = "pt", optional = true)
    public File primerPairTable = null;

    @Argument(doc="In order for this tool to design validation primers, the path to primer3 must be provided. If primer3 is in your $PATH, it will be picked up. Alternately, the environment variable PRIMER3_PATH can be set, pointing to the primer3 executable.", fullName = "primer3-path", shortName = "p3", optional = true)
    public String primer3ExePath = null;

    @Argument(doc="File to which a TSV of summary metrics should be written", fullName = "metrics-table", shortName = "mt", optional = true)
    public File metricsFile = null;

    @Argument(doc="The path to the blastn executable.  This is required for BLAST validation to be perform against putative primers. If blastn is in your $PATH, it will be picked up. Alternately, the environment variable BLASTN_PATH can be set, pointing to the blastn executable.", fullName = "blastn-path", shortName = "bn", optional = true)
    public String blastnPath = null;

    @Argument(doc="In order for this tool to use BLAST to detect validate the primers by detecting alternate binding sites, the path to a BLAST DB compiled against this reference FASTA must be provided", fullName = "blast-db-path", shortName = "bdb", optional = true)
    public String blastDatabase = null;

    @Argument(doc="If BLAST will be used, this value is passed to the -num_threads argument of blastn.", fullName = "blast-threads", shortName = "bt", optional = true)
    public Integer blastThreads = null;

    @Argument(doc="One or more files describing the insert and transgene/genome junctions. See --show-default-descriptors and --validate-descriptor-only", fullName = "insert-definition", shortName = "id", optional = true)
    public List<File> insertDescriptorFiles = null;

    @Argument(doc="The name of one of the built-in insert descriptors to use. See --write-default-descriptors for available types", fullName = "insert-name", optional = true)
    public List<String> insertNames = null;

    @Argument(doc="If provided, the tool will simply validate the insert definition files (see --insert-definition) and exit", fullName = "validate-descriptor-only", optional = true)
    public boolean validateDescriptorsOnly = false;

    @Argument(doc="If provided, the tool will write the YAML for the default descriptors to this file and exit. This can be useful as a guide for writing your own descriptors (which is likely needed).", fullName = "write-default-descriptors", optional = true)
    public File defaultDescriptorFile = null;

    @Argument(doc="If provided, this tool will inspect reads for supplemental alignments (SA tag) and parse these as well", fullName = "include-sa", shortName = "sa", optional = true)
    public boolean includeSupplementalAlignments = false;

    @Override
    public boolean requiresReference() {
        return !validateDescriptorsOnly;
    }

    private List<InsertDescriptor> INSERT_DESCRIPTORS = new ArrayList<>();

    private static final ObjectMapper OM = new ObjectMapper(new YAMLFactory());

    private InsertDescriptor parseInsertType(File yml) {
        try {
            return OM.readValue(yml, InsertDescriptor.class);
        }
        catch (Exception e) {
            throw new UserException.BadInput("Unable to parse file: " + yml.getName() + ", error was: " + e.getMessage());
        }
    }

    private enum DefaultDescriptorType {
        lentivirus("Lentivirus.yml"),
        piggybac("pENTR-PB511-Repro.yml");

        private String fileName;

        DefaultDescriptorType(String fileName) {
            this.fileName = fileName;
        }

        public InputStream getStream() {
            Resource yml = new Resource(fileName, TagPcrSummary.class);

            return yml.getResourceContentsAsStream();
        }

    }
    private InsertDescriptor getDefaultType(String name) {
        try {
            DefaultDescriptorType d = DefaultDescriptorType.valueOf(name.toLowerCase());
            try {
                return OM.readValue(d.getStream(), InsertDescriptor.class);
            }
            catch (Exception e) {
                throw new UserException.BadInput("Unable to parse insert definition: " + name + ", error was: " + e.getMessage());
            }
        }
        catch (IllegalArgumentException e) {
            throw new UserException.BadInput("Unknown insert type: " + name);
        }
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        if (defaultDescriptorFile != null) {
            IOUtil.assertFileIsWritable(defaultDescriptorFile);

            try (PrintStream writer = new PrintStream(IOUtil.openFileForWriting(defaultDescriptorFile))) {
                writer.println("Below are the built-in insert descriptors. These can serve as examples for creating a custom descriptor. Multiple examples are below, but the tool will expect one per file if you create you own.");
                writer.println("");
                for (DefaultDescriptorType d : DefaultDescriptorType.values()) {
                    writer.println("Example: " + d.name());
                    writer.println("");

                    IOUtil.copyStream(d.getStream(), writer);
                    writer.println("");
                    writer.println("");
                }
            }

            logger.info("Default insert descriptors written to: " + defaultDescriptorFile.getPath());
            return;
        }

        for (File insertType : insertDescriptorFiles) {
            INSERT_DESCRIPTORS.add(parseInsertType(insertType));
        }

        for (String name : insertNames) {
            INSERT_DESCRIPTORS.add(getDefaultType(name));
        }

        // The block above will error if not valid
        if (validateDescriptorsOnly) {
            logger.info("The insert descriptors were valid");
            return;
        }
        else {
            if (inputBam == null) {
                throw new UserException.CouldNotReadInputFile("BAM not provided");
            }

            IOUtil.assertFileIsReadable(inputBam);

            if (outputTsv == null) {
                throw new UserException.BadInput("The argument --output-table is required.");
            }

            IOUtil.assertFileIsWritable(outputTsv);
        }

        if (blastDatabase != null) {
            File blastTest = new File(blastDatabase + ".nhr");
            if (!blastTest.exists()) {
                throw new UserException.BadInput("Invalid BLAST DB path.  Expected this location to contain many files, including: " + blastTest.getPath());
            }

            blastnPath = setExePath(blastnPath, "blastn", "BLASTN_PATH");
        }
        else {
            logger.info("No BLAST database provided, will not validate primers");
        }

        primer3ExePath = setExePath(primer3ExePath, "primer3-core", "PRIMER3_PATH");
        if (doGenerateAmplicons() && primer3ExePath == null) {
            throw new UserException.BadInput("Primer generation was implied, but primer3 was not found");
        }
    }

    private String setExePath(String initialValue, String exe, String environmentVariable) {
        if (SystemUtils.IS_OS_WINDOWS) {
            exe = exe + ".exe";
        }

        // Preferentially take user-supplied argument
        if (initialValue != null) {
            if (!new File(initialValue).exists()) {
                throw new UserException.BadInput("File not found: " + initialValue);
            }

            return initialValue;
        }

        // Next try environment:
        String env = StringUtils.trimToNull(System.getenv(environmentVariable));
        if (env != null && new File(env).exists()) {
            return env;
        }

        // Finally iterate PATH
        return findExecutableOnPath(exe);
    }

    private String findExecutableOnPath(String name) {
        for (String dirname : System.getenv("PATH").split(File.pathSeparator)) {
            File file = new File(dirname, name);
            if (file.isFile() && file.canExecute()) {
                return file.getAbsolutePath();
            }
        }

        return null;
    }

    private boolean exitWithoutRunning() {
        return validateDescriptorsOnly || defaultDescriptorFile != null;
    }

    private boolean doGenerateAmplicons() {
        return (blastDatabase != null && primerPairTable != null) || outputGenbank != null;
    }

    private Map<String, Integer> primaryAlignmentsMatchingInsert = new HashMap<>();
    private Map<String, Integer> secondaryAlignmentsMatchingInsert = new HashMap<>();

    private void inspectForInsert(SAMRecord rec) {
        final String read = rec.getReadString().toUpperCase();
        INSERT_DESCRIPTORS.forEach(id -> {
            for (String ss : id.getAllBackboneSearchStrings()) {
                Integer isMatch = SequenceMatcher.fuzzyMatch(ss.toUpperCase(), read, id.getBackboneSearchEditDistance());
                if (isMatch != null) {
                    if (rec.isSecondaryOrSupplementary()) {
                        int total = secondaryAlignmentsMatchingInsert.getOrDefault(id.getName(), 0);
                        total++;
                        secondaryAlignmentsMatchingInsert.put(id.getName(), total);
                    }
                    else {
                        int total = primaryAlignmentsMatchingInsert.getOrDefault(id.getName(), 0);
                        total++;
                        primaryAlignmentsMatchingInsert.put(id.getName(), total);
                    }

                    break;
                }
            }
        });
    }

    @Override
    public void traverse() {
        if (exitWithoutRunning()) {
            return;
        }

        NumberFormat pf = NumberFormat.getPercentInstance();
        pf.setMaximumFractionDigits(6);

        SamReaderFactory fact = SamReaderFactory.makeDefault();
        fact.validationStringency(ValidationStringency.SILENT);
        fact.referenceSequence(referenceArguments.getReferencePath());

        int totalAlignments = 0;
        int uniqueReads = 0;
        int numReadsSpanningJunction = 0;
        int splitAlignments = 0;
        File bam = ensureQuerySorted(inputBam);

        Map<String, JunctionMatch> totalMatches = new HashMap<>();
        Map<String, Number> metricsMap = new HashMap<>();

        int reverseReadsSkipped = 0;

        try (SamReader bamReader = fact.open(bam); SAMRecordIterator it = bamReader.iterator()) {

            List<SAMRecord> alignmentsForRead = new ArrayList<>();
            while (it.hasNext()) {
                SAMRecord rec = it.next();
                if (rec.getReadUnmappedFlag()) {
                    inspectForInsert(rec);
                }

                // Skip reverse reads
                if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                    reverseReadsSkipped++;
                    continue;
                }

                //If this read doesnt match the prior set, process these and clear alignmentsForRead
                if (!alignmentsForRead.isEmpty() && !alignmentsForRead.get(0).getReadName().equals(rec.getReadName())){
                    uniqueReads++;
                    numReadsSpanningJunction += processAlignmentsForRead(alignmentsForRead, totalMatches);
                }

                totalAlignments++;
                alignmentsForRead.add(rec);

                if (includeSupplementalAlignments && rec.hasAttribute("SA")) {
                    String sa = StringUtils.trimToNull(rec.getStringAttribute("SA"));
                    if (sa != null)
                    {
                        for (String alignment : sa.split(";"))
                        {
                            String[] parts = alignment.split(",");
                            SAMRecord newRec = rec.deepCopy();
                            newRec.setReferenceIndex(rec.getHeader().getSequenceIndex(parts[0]));
                            newRec.setAlignmentStart(Integer.parseInt(parts[1]));
                            newRec.setReadNegativeStrandFlag("-".equals(parts[2]));
                            if (rec.getReadNegativeStrandFlag() != newRec.getReadNegativeStrandFlag()) {
                                newRec.reverseComplement();
                            }
                            newRec.setCigar(TextCigarCodec.decode(parts[3]));
                            newRec.setMappingQuality(Integer.parseInt(parts[4]));

                            alignmentsForRead.add(newRec);
                            splitAlignments++;
                        }
                    }
                }
            }

            //ensure we capture final read
            if (!alignmentsForRead.isEmpty()) {
                uniqueReads++;
                numReadsSpanningJunction += processAlignmentsForRead(alignmentsForRead, totalMatches);

            }
        }
        catch (IOException e)
        {
            throw new GATKException(e.getMessage(), e);
        }

        double pct = uniqueReads == 0 ? 0 : numReadsSpanningJunction / (double)uniqueReads;
        logger.info("Total reads spanning a junction: " + numReadsSpanningJunction + " of " + uniqueReads + " (" + pf.format(pct) + ")");

        metricsMap.put("TotalAlignments", totalAlignments);
        metricsMap.put("SplitAlignments", splitAlignments);
        metricsMap.put("DistinctReads", uniqueReads);
        metricsMap.put("NumReadsSpanningJunction", numReadsSpanningJunction);
        metricsMap.put("PctReadsSpanningJunction", pct);

        Set<String> inserts = new TreeSet<>();
        inserts.addAll(primaryAlignmentsMatchingInsert.keySet());
        inserts.addAll(secondaryAlignmentsMatchingInsert.keySet());
        logger.info("Total primary/secondary alignments matching insert/transgene:");
        int totalPrimaryAlignmentsMatchingInsert = 0;
        int totalSecondaryAlignmentsMatchingInsert = 0;
        for (String key : inserts) {
            logger.info(key + ", primary: " + primaryAlignmentsMatchingInsert.getOrDefault(key, 0));
            totalPrimaryAlignmentsMatchingInsert += primaryAlignmentsMatchingInsert.getOrDefault(key, 0);

            logger.info(key + ", secondary: " + secondaryAlignmentsMatchingInsert.getOrDefault(key, 0));
            totalSecondaryAlignmentsMatchingInsert += secondaryAlignmentsMatchingInsert.getOrDefault(key, 0);
        }

        metricsMap.put("TotalPrimaryAlignmentsMatchingInsert", totalPrimaryAlignmentsMatchingInsert);
        metricsMap.put("FractionPrimaryAlignmentsMatchingInsert", (double)totalPrimaryAlignmentsMatchingInsert / uniqueReads);
        metricsMap.put("TotalSecondaryAlignmentsMatchingInsert", totalSecondaryAlignmentsMatchingInsert);

        int totalMinusStrand = 0;
        List<DNASequence> amplicons = new ArrayList<>();

        if (!totalMatches.isEmpty()) {
            NumberFormat format = DecimalFormat.getNumberInstance();
            format.setMaximumFractionDigits(6);
            try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(outputTsv), '\t', CSVWriter.NO_QUOTE_CHARACTER); CSVWriter primerWriter = primerPairTable == null ? null : new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(primerPairTable), '\t', CSVWriter.NO_QUOTE_CHARACTER); ReferenceSequenceFile refSeq = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceArguments.getReferencePath())) {
                writer.writeNext(new String[]{"SiteName", "JunctionName", "Orientation", "Chr", "Position", "Strand", "Total", "Fraction"});
                if (primerWriter != null) {
                    primerWriter.writeNext(new String[]{"SiteName", "JunctionName", "Orientation", "Primer-F-Name", "Primer-F-Seq", "Primer-F-Start", "Primer-R-Name", "Primer-R-Seq", "Primer-R-Start"});
                }
                Map<String, ReferenceSequence> refMap = new HashMap<>();
                List<String> sortedKeys = new ArrayList<>(totalMatches.keySet());
                sortedKeys.sort(ComparatorUtils.naturalComparator());

                AtomicInteger totalHits = new AtomicInteger();
                totalMatches.forEach((x, y) -> totalHits.getAndAdd(y.totalReads));

                int totalOutput = 0;
                for (String key : sortedKeys) {
                    JunctionMatch jm = totalMatches.get(key);

                    if (jm.totalReads < MIN_ALIGNMENTS) {
                        continue;
                    }

                    double fraction = (double)jm.totalReads / totalHits.get();
                    if (fraction < MIN_FRACTION) {
                        continue;
                    }

                    totalOutput++;

                    ReferenceSequence ref = refMap.get(jm.contigName);
                    if (ref == null) {
                        ref = refSeq.getSequence(jm.contigName);
                        refMap.put(jm.contigName, ref);
                    }

                    if (jm.isTransgeneInverted()) {
                        totalMinusStrand++;
                    }

                    String siteName = "Site" + totalOutput;
                    writer.writeNext(new String[]{siteName, jm.jd.getName(), jm.getTransgeneOrientation(), jm.contigName, String.valueOf(jm.matchStart), jm.getTransgeneStrand(), String.valueOf(jm.totalReads), format.format(fraction)});

                    if (READS_PER_SITE > 0 && !jm.representativeRecords.isEmpty()) {
                        File fastaOut = new File(outputTsv.getParentFile(), siteName + "." + jm.jd.getName() + "." + jm.getTransgeneOrientation() + ".fasta.gz");
                        try (BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openGzipFileForWriting(fastaOut.toPath())))) {
                            for (SAMRecord r : jm.representativeRecords) {
                                out.write(">" + r.getReadName() + (r.getReadNegativeStrandFlag() ? "-R" : "") + "\n");
                                String readString = StringUtil.bytesToString(r.getReadBases());
                                if (r.getReadNegativeStrandFlag())
                                    readString = SequenceUtil.reverseComplement(readString);

                                out.write(readString + "\n");
                            }
                        }
                    }

                    if (doGenerateAmplicons()) {
                        Pair<DNASequence, DNASequence> ampliconPair = jm.getAmplicons(refSeq, refMap, siteName, primerWriter, outputTsv.getParentFile(), primer3ExePath);
                        amplicons.add(ampliconPair.getKey());
                        amplicons.add(ampliconPair.getValue());
                    }
                }

                logger.info("Total junctions output: " + totalOutput);
                logger.info("Plus Strand: " + (totalOutput - totalMinusStrand) + " / Minus: " + totalMinusStrand);

                metricsMap.put("TotalIntegrationSitesOutput", totalOutput);
                metricsMap.put("IntegrationSitesOutputMinusStrand", totalMinusStrand);
                metricsMap.put("IntegrationSitesOutputPlusStrand", (totalOutput - totalMinusStrand));
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            if (blastDatabase != null && primerPairTable != null) {
                runBlastN(primerPairTable, amplicons);
            }
        }
        else {
            logger.info("there were no passing alignments");
        }

        //Now make genbank output:
        if (outputGenbank != null && !amplicons.isEmpty()) {
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

    private int processAlignmentsForRead(List<SAMRecord> alignmentsForRead, Map<String, JunctionMatch> totalMatches) {
        Map<String, JunctionMatch> matches = new HashMap<>();
        alignmentsForRead.forEach(rec -> {
            if (rec.isSecondaryOrSupplementary() && !includeSupplementalAlignments) {
                return;
            }

            for (InsertDescriptor id : INSERT_DESCRIPTORS) {
                for (InsertJunctionDescriptor jd : id.getJunctions()) {
                    Map<String, JunctionMatch> hits = jd.getMatches(rec, id, logger, READS_PER_SITE);
                    matches.putAll(hits);
                }
            }

            inspectForInsert(rec);
        });

        alignmentsForRead.clear();
        if (!matches.isEmpty()) {
            matches.forEach((x, y) -> {
                if (totalMatches.containsKey(x)) {
                    totalMatches.get(x).addRead(y.representativeRecords);
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
            logger.info("writing to file: " + querySorted.getPath());

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

    public enum END_TYPE {
        downstream(),
        upstream()
    }

    protected static class JunctionMatch {
        private final Logger log;
        private final InsertDescriptor insertDescriptor;
        private final InsertJunctionDescriptor jd;

        private final String contigName;
        private final int matchStart;
        private final boolean hitOrientationReversed;
        private int totalReads = 1;
        private List<SAMRecord> representativeRecords = new ArrayList<>();
        private final int maxRecordsToStore;

        public JunctionMatch(Logger log, InsertDescriptor insertDescriptor, InsertJunctionDescriptor jd, String contigName, int matchStart, boolean hitOrientationReversed, SAMRecord rec, int maxRecordsToStore) {
            this.log = log;
            this.insertDescriptor = insertDescriptor;
            this.jd = jd;
            this.contigName = contigName;
            this.matchStart = matchStart;
            this.hitOrientationReversed = hitOrientationReversed;
            this.maxRecordsToStore = maxRecordsToStore;

            if (maxRecordsToStore > 0 && rec != null)
                representativeRecords.add(rec);
        }

        public String getKey() {
            return jd.getName() + "<>" + contigName + "<>" + matchStart + "<>" + hitOrientationReversed;
        }

        public void addRead(List<SAMRecord> recs) {
            totalReads += 1;
            if (representativeRecords.size() < maxRecordsToStore && recs != null) {
                representativeRecords.addAll(recs);
            }
        }

        public boolean isTransgeneInverted()
        {
            if (hitOrientationReversed) {
                return !jd.isInvertHitOrientation();
            }
            else {
                return jd.isInvertHitOrientation();
            }
        }

        public String getTransgeneStrand()
        {
            return isTransgeneInverted() ? "-" : "+";
        }

        public String getTransgeneOrientation()
        {
            return isTransgeneInverted() ? "Minus" : "Plus";
        }

        public Pair<DNASequence, DNASequence> getAmplicons(ReferenceSequenceFile refSeq, Map<String, ReferenceSequence> refMap, String siteName, CSVWriter primerWriter, File outDir, String primer3ExePath) {
            ReferenceSequence ref = refMap.get(contigName);
            if (ref == null) {
                ref = refSeq.getSequence(contigName);
                refMap.put(contigName, ref);
            }

            Interval afterInterval = new Interval(ref.getName(), matchStart, Math.min(matchStart + 1000, ref.length()));
            String genomeAfter = ref.getBaseString().substring(afterInterval.getStart() - 1, afterInterval.getEnd());

            Interval beforeInterval = new Interval(ref.getName(), Math.max(1, matchStart - 1000), matchStart);
            String genomeBefore = ref.getBaseString().substring(beforeInterval.getStart() - 1, beforeInterval.getEnd());

            String insertUpstreamRegion = "";
            if (insertDescriptor.getInsertUpstreamRegion() != null) {
                insertUpstreamRegion = insertDescriptor.getInsertUpstreamRegion().getSequence();
            }

            String insertDownstreamRegion = "";
            if (insertDescriptor.getInsertDownstreamRegion() != null) {
                insertDownstreamRegion = insertDescriptor.getInsertDownstreamRegion().getSequence();
            }

            try {
                String id = siteName;
                String idFull = contigName + ":" + matchStart + ":" + getTransgeneOrientation();

                //Note: this is a limitation in biojava
                if (id.length() > 15) {
                    log.error("A site name greater than 15 characters was passed: " + id + ".  This will be truncated to comply with genbank format");
                    id = id.substring(0, 15);
                }

                String orientationSuffix = isTransgeneInverted() ? "m" : "";

                String effectiveRegionUpstream = isTransgeneInverted() ? SequenceUtil.reverseComplement(insertDownstreamRegion) : insertUpstreamRegion;
                String effectiveRegionUpstreamLabel = isTransgeneInverted() ? insertDescriptor.getFeatureLabel(END_TYPE.downstream) + "-RC" : insertDescriptor.getFeatureLabel(END_TYPE.upstream);

                String effectiveRegionDownstream = isTransgeneInverted() ? SequenceUtil.reverseComplement(insertUpstreamRegion) : insertDownstreamRegion;
                String effectiveRegionDownstreamLabel = isTransgeneInverted() ? insertDescriptor.getFeatureLabel(END_TYPE.upstream) + "-RC" : insertDescriptor.getFeatureLabel(END_TYPE.downstream);

                DNASequence genomicJunctionUpstream = new DNASequence(genomeBefore + effectiveRegionUpstream);
                String ampliconName = id + "-" + effectiveRegionUpstreamLabel;
                String accession = ampliconName + orientationSuffix;
                if (accession.length() > 15) {
                    //a limitation of BioJava
                    accession = accession.substring(0, 15);
                }
                genomicJunctionUpstream.setAccession(new AccessionID(accession));
                genomicJunctionUpstream.setDescription(idFull);
                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(effectiveRegionUpstreamLabel, "Vector", jd.getName(), "");
                tf.setLocation(new SequenceLocation<>(genomeBefore.length() + 1, genomicJunctionUpstream.getBioEnd(), genomicJunctionUpstream, isTransgeneInverted() ? Strand.NEGATIVE : Strand.POSITIVE));
                genomicJunctionUpstream.addFeature(tf);
                runPrimer3(siteName, effectiveRegionUpstreamLabel, genomicJunctionUpstream, genomeBefore.length(), outDir, primerWriter, primer3ExePath);

                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tfG = new TextFeature<>(beforeInterval.getContig() + ":" + beforeInterval.getStart(), "Genome", beforeInterval.toString(), "");
                tfG.setLocation(new SequenceLocation<>(1, genomeBefore.length(), genomicJunctionUpstream, Strand.POSITIVE));
                genomicJunctionUpstream.addFeature(tfG);

                DNASequence seqAfterInsert = new DNASequence(effectiveRegionDownstream + genomeAfter);
                String ampliconName2 = id + "-" + effectiveRegionDownstreamLabel;
                String accession2 = ampliconName2 + orientationSuffix;
                if (accession2.length() > 15) {
                    //a limitation of BioJava
                    accession2 = accession2.substring(0, 15);
                }
                seqAfterInsert.setAccession(new AccessionID(accession2));
                seqAfterInsert.setDescription(idFull);
                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf2 = new TextFeature<>(effectiveRegionDownstreamLabel, "Vector", jd.getName(), "");
                tf2.setLocation(new SequenceLocation<>(1, effectiveRegionDownstream.length(), seqAfterInsert, isTransgeneInverted() ? Strand.NEGATIVE : Strand.POSITIVE));
                seqAfterInsert.addFeature(tf2);
                runPrimer3(siteName, effectiveRegionDownstreamLabel, seqAfterInsert, effectiveRegionDownstream.length(), outDir, primerWriter, primer3ExePath);

                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tfG2 = new TextFeature<>(afterInterval.getContig() + ":" + afterInterval.getStart(), "Genome", afterInterval.toString(), afterInterval.toString());
                tfG2.setLocation(new SequenceLocation<>(effectiveRegionDownstream.length() + 1, seqAfterInsert.getBioEnd(), seqAfterInsert, Strand.POSITIVE));
                seqAfterInsert.addFeature(tfG2);

                //Add known primers as features:
                if (!insertDescriptor.getInternalPrimers().isEmpty()) {
                    for (SequenceDescriptor p : insertDescriptor.getInternalPrimers()) {
                        addPrimerIfPresent(p, genomicJunctionUpstream);
                        addPrimerIfPresent(p, seqAfterInsert);
                    }
                }

                return Pair.of(genomicJunctionUpstream, seqAfterInsert);

            } catch (CompoundNotFoundException e) {
                throw new GATKException(e.getMessage(), e);
            }
        }

        private String getTruncatedLabel(String val, boolean addRCSuffix) {
            int max = addRCSuffix ? 12 : 15;

            if (val.length() > max) {
                val = val.substring(0, max);
            }

            return val + (addRCSuffix ? "-RC" : "");
        }

        private void addPrimerIfPresent(SequenceDescriptor primer, DNASequence seq) {
            int pos = seq.getSequenceAsString().indexOf(primer.getSequence());
            if (pos > -1) {
                pos++; //1-based
                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(getTruncatedLabel(primer.getName(), false), "Vector Primer", primer.getName(), "Existing Primer: " + primer.getName());
                tf.setLocation(new SequenceLocation<>(pos, pos + primer.getSequence().length(), seq, Strand.POSITIVE));
                seq.addFeature(tf);
            }

            //Also test RC:
            try {
                String primerSeqRC = new DNASequence(primer.getSequence()).getReverseComplement().getSequenceAsString();
                int posRC = seq.getSequenceAsString().indexOf(primerSeqRC);
                if (posRC > -1) {
                    posRC++; //1-based
                    TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(getTruncatedLabel(primer.getName(), true), "Vector Primer", getTruncatedLabel(primer.getName(), true), "Existing Primer: " + primer.getName() + "-RC");
                    tf.setLocation(new SequenceLocation<>(posRC, posRC + primer.getSequence().length(), seq, Strand.NEGATIVE));
                    seq.addFeature(tf);
                }
            } catch (CompoundNotFoundException e) {
                throw new GATKException(e.getMessage(), e);
            }
        }

        private void runPrimer3(String siteName, String junctionName, DNASequence sequence, int junctionSite, File outDir, CSVWriter primerWriter, String primer3ExePath) {
            log.info("Running primer3 for: " + siteName + ", " + junctionName);

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
                                "PRIMER_PRODUCT_SIZE_RANGE=250-600 601-1000 " + "\n" +
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
                log.info("primer3 had a non-zero exit");
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

            if (outputMap.containsKey("PRIMER_ERROR")) {
                log.info("primer3 error: " + outputMap.get("PRIMER_ERROR"));
            }

            if (!outputMap.containsKey("PRIMER_PAIR_NUM_RETURNED")) {
                log.info("primer3 returned no pairs");
                return;
            }

            output.delete();

            int numPairs = Integer.parseInt(outputMap.get("PRIMER_PAIR_NUM_RETURNED"));
            log.info("total candidate pairs: " + numPairs);

            for (int i = 0; i < numPairs; i++) {
                String primer1Seq = outputMap.get("PRIMER_LEFT_" + i + "_SEQUENCE");
                int primer1Loc = Integer.parseInt(outputMap.get("PRIMER_LEFT_" + i).split(",")[0]) + 1;

                String primer2Seq = outputMap.get("PRIMER_RIGHT_" + i + "_SEQUENCE");
                int primer2Loc = Integer.parseInt(outputMap.get("PRIMER_RIGHT_" + i).split(",")[0]) + 1;

                String pairNameShort = junctionName + "-" + (i + 1);
                String pairName = siteName + "-" + pairNameShort;
                if (primerWriter != null) {
                    primerWriter.writeNext(new String[]{siteName, junctionName, getTransgeneOrientation(), pairName + "F", primer1Seq, String.valueOf(primer1Loc), pairName + "R", primer2Seq, String.valueOf(primer2Loc)});
                }

                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf = new TextFeature<>(pairNameShort + "F", "Primer3", pairName + "F", pairName + "F");
                tf.setLocation(new SequenceLocation<>(primer1Loc, primer1Loc + primer1Seq.length() - 1, sequence, Strand.POSITIVE));
                sequence.addFeature(tf);

                TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> tf2 = new TextFeature<>(pairNameShort + "R", "Primer3", pairName + "R", pairName + "R");
                tf2.setLocation(new SequenceLocation<>(primer2Loc - primer2Seq.length() + 1, primer2Loc, sequence, Strand.NEGATIVE));
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

            logger.info("Duplicate primers collapsed before BLAST: " + primersSkipped + " of " + (idx *2) + ".  Unique primers to BLAST: " + idxToPrimer.size());
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
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        blastInput.delete();
        blastOutput.delete();

        //Now write updated output:
        Set<String> primerNamesToFlag = new HashSet<>();
        try (CSVWriter primerWriter = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(primer3Table), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            primerWriter.writeNext(new String[]{"SiteName", "JunctionName", "Orientation", "Primer-F-Name", "Primer-F-Seq", "Primer-F-Start", "Primer-R-Name", "Primer-R-Seq", "Primer-R-Start", "ForwardPassedBlast", "ReversePassedBlast", "BothPassedBlast"});
            for (String[] line : primerLines) {
                boolean forwardPass = !failedPrimers.containsKey(line[4].toUpperCase());
                boolean reversePass = !failedPrimers.containsKey(line[7].toUpperCase());
                List<String> toWrite = new ArrayList<>(Arrays.asList(line));
                toWrite.add(forwardPass ? "Y" : "");
                toWrite.add(reversePass ? "Y" : "");
                toWrite.add((forwardPass && reversePass) ? "Y" : "");
                primerWriter.writeNext(toWrite.toArray(new String[toWrite.size()]));

                if (!forwardPass) {
                    primerNamesToFlag.add(line[3]);
                }
                
                if (!reversePass) {
                    primerNamesToFlag.add(line[6]);
                }
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        logger.info("Total primer paired discarded due to multiple BLAST hits: " + (primerNamesToFlag.size() / 2) + " of " + primerLines.size());

        //prune failed from genbank output:
        if (!primerNamesToFlag.isEmpty()) {
            amplicons.listIterator().forEachRemaining(seq -> {
                List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = new ArrayList<>(seq.getFeatures());
                features.listIterator().forEachRemaining(feature -> {
                    if (primerNamesToFlag.contains(feature.getDescription())) {
                        feature.setType("Fail:" + feature.getType());
                    }
                });
            });
        }
    }
}


