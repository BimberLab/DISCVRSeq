package com.github.discvrseq.tools;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.util.SequenceMatcher;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.collections4.ComparatorUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.TextFeature;
import org.biojava.nbio.core.sequence.io.GenbankWriter;
import org.biojava.nbio.core.sequence.io.GenericGenbankHeaderFormat;
import org.biojava.nbio.core.sequence.location.SequenceLocation;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.*;
import java.text.NumberFormat;
import java.util.*;

/**
 * This tool was originally created as part of an annotation pipeline for non-human data.  The input VCF from another species (or genome build) would be lifted to the human genome and annotated
 * in those coordinates. Lifeover must be performed using Picard Tools LiftoverVcf, which annotates lifted variants with the values for ORIGNAL_CONTIG, ORGINAL_START and ORIGINAL_ALLELE.  This tool
 * reads one of these lifted VCFs, and writes a new sorted VCF in which the original coordinates are restored.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar TagPcrSummary \
 *     -R currentGenome.fasta \
 *     -I myBam.bam \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = DiscvrSeqInternalProgramGroup.class
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
    }

    @Override
    public void traverse() {
        NumberFormat pf = NumberFormat.getPercentInstance();
        pf.setMaximumFractionDigits(6);

        SamReaderFactory fact = SamReaderFactory.makeDefault();
        fact.validationStringency(ValidationStringency.SILENT);
        fact.referenceSequence(referenceArguments.getReferencePath());

        int totalReads = 0;
        int numReadsSpanningJunction = 0;
        File bam = ensureQuerySorted(inputBam);

        Map<String, JunctionMatch> totalMatches = new HashMap<>();

        try (SamReader bamReader = fact.open(bam); SAMRecordIterator it = bamReader.iterator()) {

            List<SAMRecord> alignmentsForRead = new ArrayList<>();
            while (it.hasNext()) {
                SAMRecord rec = it.next();
                if (rec.getReadUnmappedFlag()) {
                    continue;
                }

                // Skip reverse reads
                if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                    continue;
                }

                if (alignmentsForRead.isEmpty() || alignmentsForRead.get(0).getReadName().equals(rec.getReadName())) {
                    alignmentsForRead.add(rec);
                } else {
                    totalReads++;
                    numReadsSpanningJunction += processAlignmentsForRead(alignmentsForRead, totalMatches);
                }
            }

            //ensure we capture final read
            if (!alignmentsForRead.isEmpty()) {
                totalReads++;
                numReadsSpanningJunction += processAlignmentsForRead(alignmentsForRead, totalMatches);

            }
        }
        catch (IOException e)
        {
            throw new GATKException(e.getMessage(), e);
        }

        double pct = totalReads == 0 ? 0 : numReadsSpanningJunction / (double)totalReads;
        logger.info("Total reads spanning a junction: " + numReadsSpanningJunction + " of " + totalReads + " (" + pf.format(pct) + ")");

        int totalReverse = 0;
        if (!totalMatches.isEmpty()) {
            try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(outputTsv), '\t', CSVWriter.NO_QUOTE_CHARACTER); CSVWriter primerWriter = primerPairTable == null ? null : new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(primerPairTable), '\t', CSVWriter.NO_QUOTE_CHARACTER);ReferenceSequenceFile refSeq = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceArguments.getReferencePath())) {
                writer.writeNext(new String[]{"SiteName", "JunctionName", "Chr", "Position", "Strand", "Total"});
                if (primerWriter != null) {
                    primerWriter.writeNext(new String[]{"SiteName", "JunctionName", "Primer-F", "Primer-F-Start", "Primer-R", "Primer-R-Start"});
                }
                Map<String, ReferenceSequence> refMap = new HashMap<>();
                List<String> sortedKeys = new ArrayList<>(totalMatches.keySet());
                sortedKeys.sort(ComparatorUtils.naturalComparator());

                int totalOutput = 0;

                List<DNASequence> amplicons = new ArrayList<>();

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
                        writer.writeNext(new String[]{siteName, jm.jd.junctionName, jm.contigName, String.valueOf(jm.alignmentStart), (jm.matchIsNegativeStrand ? "-" : "+"), String.valueOf(jm.totalReads)});

                        Pair<DNASequence, DNASequence> ampliconPair = jm.getAmplicons(refSeq, refMap, siteName, primerWriter);
                        amplicons.add(ampliconPair.getKey());
                        amplicons.add(ampliconPair.getValue());

                    } else {
                        logger.info("Skipping position: " + key + ", " + jm.totalReads + " reads");
                    }
                }

                logger.info("Total junctions output: " + totalOutput);
                logger.info("Forward / Reverse: " + (totalOutput - totalReverse) + " / " + totalReverse);

                //Now make genbank output:
                if (outputGenbank != null) {
                    try (BufferedOutputStream os = new BufferedOutputStream(IOUtil.openFileForWriting(outputGenbank))) {
                        GenbankWriter<DNASequence, NucleotideCompound> genbankWriter = new GenbankWriter<>(os, amplicons, new GenericGenbankHeaderFormat<>("LINEAR"));
                        genbankWriter.process();
                    }
                }
            }
            catch (Exception e)
            {
                throw new GATKException(e.getMessage(), e);
            }
        }
    }

    private int processAlignmentsForRead(List<SAMRecord> alignmentsForRead, Map<String, JunctionMatch> totalMatches) {
        Map<String, JunctionMatch> matches = new HashMap<>();
        alignmentsForRead.forEach(rec -> {
            if (rec.isSecondaryOrSupplementary()) {
                logger.info(rec.getReferenceName());
            }
            else {
                for (InsertDescriptor id : INSERT_DESCRIPTORS) {
                    for (InsertJunctionDescriptor jd : id.junctions) {
                        matches.putAll(jd.getMatches(rec));
                    }
                }
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

        public Map<String, JunctionMatch> getMatches(SAMRecord rec) {
            Map<String, JunctionMatch> matches = new HashMap<>();

            //Forward orientation.
            for (String query : getSearchStrings(false)) {
                Integer match = SequenceMatcher.fuzzyMatch(query, rec.getReadString().toUpperCase(), mismatchesAllowed);
                if (match != null) {
                    //let the read deal with clipping
                    JunctionMatch jm = new JunctionMatch(this, rec.getContig(), rec.getAlignmentStart(), false);
                    matches.put(jm.getKey(), jm);
                }
            }

            //Reverse:
            for (String query : getSearchStrings(rec.getReadNegativeStrandFlag())) {
                Integer match = SequenceMatcher.fuzzyMatch(query, rec.getReadString().toUpperCase(), mismatchesAllowed);
                if (match != null) {
                    JunctionMatch jm = new JunctionMatch(this, rec.getContig(), rec.getAlignmentEnd(), true);
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
        private int totalReads = 1;

        public JunctionMatch(InsertJunctionDescriptor jd, String contigName, int alignmentStart, boolean matchIsNegativeStrand) {
            this.jd = jd;
            this.contigName = contigName;
            this.alignmentStart = alignmentStart;
            this.matchIsNegativeStrand = matchIsNegativeStrand;
        }

        public String getKey() {
            return jd.junctionName + "<>" + contigName + "<>" + alignmentStart + "<>" + matchIsNegativeStrand;
        }

        public void addRead() {
            totalReads += 1;
        }

        public Pair<DNASequence, DNASequence> getAmplicons(ReferenceSequenceFile refSeq, Map<String, ReferenceSequence> refMap, String siteName, CSVWriter primerWriter ) {
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
                    refMap.put(insertRef.getName(), ref);
                }

                insertRegionProximal = ref.getBaseString().substring(jd.insertRegionProximal.getStart() - 1, jd.insertRegionProximal.getEnd());
            }

            String insertRegionDistal = "";
            if (jd.insertRegionDistal != null) {
                ReferenceSequence insertRef = refMap.get(jd.insertRegionDistal.getContig());
                if (insertRef == null) {
                    insertRef = refSeq.getSequence(jd.insertRegionDistal.getContig());
                    refMap.put(insertRef.getName(), ref);
                }

                insertRegionDistal = ref.getBaseString().substring(jd.insertRegionDistal.getStart() - 1, jd.insertRegionDistal.getEnd());
            }

            DNASequence seq1;
            DNASequence seq2;
            try {
                String id = siteName;
                String idFull = contigName + ":" + alignmentStart + ":" + (matchIsNegativeStrand ? "Minus" : "Plus");

                //TODO: fix limitation
                if (id.length() > 16) {
                    id = id.substring(0, 16);
                }

                if (matchIsNegativeStrand) {
                    seq1 = new DNASequence(SequenceUtil.reverseComplement(insertRegionDistal) + after);
                    String ampliconName = id + "-" + jd.getFeatureLabel(END_TYPE.distal);
                    seq1.setAccession(new AccessionID(ampliconName));
                    seq1.setDescription(idFull);

                    TextFeature tf = new TextFeature<>(jd.getFeatureLabel(END_TYPE.distal) + "-RC", "Vector", jd.junctionName, jd.insertRegionDistal.toString());
                    tf.setLocation(new SequenceLocation<>(1, insertRegionDistal.length(), seq1, Strand.NEGATIVE));
                    seq1.addFeature(tf);
                    runPrimer3(ampliconName, jd.getFeatureLabel(END_TYPE.distal), seq1, insertRegionDistal.length(), outputTsv.getParentFile(), primerWriter);

                    TextFeature tfG = new TextFeature<>(afterInterval.toString(), "Genome", afterInterval.toString(), afterInterval.toString());
                    tfG.setLocation(new SequenceLocation<>(insertRegionDistal.length() + 1, seq1.getBioEnd(), seq1, Strand.NEGATIVE));
                    seq1.addFeature(tfG);

                    seq2 = new DNASequence(before + SequenceUtil.reverseComplement(insertRegionProximal));
                    String ampliconName2 = id + "-" + jd.getFeatureLabel(END_TYPE.promimal);
                    seq2.setAccession(new AccessionID(ampliconName2));
                    seq2.setDescription(idFull);
                    TextFeature tf2 = new TextFeature<>(jd.getFeatureLabel(END_TYPE.promimal) + "-RC", "Vector", jd.junctionName, jd.insertRegionProximal.toString());
                    tf2.setLocation(new SequenceLocation<>(before.length() + 1, seq2.getBioEnd(), seq2, Strand.NEGATIVE));
                    seq2.addFeature(tf2);
                    runPrimer3(ampliconName2, jd.getFeatureLabel(END_TYPE.promimal), seq2, before.length(), outputTsv.getParentFile(), primerWriter);
                } else {
                    seq1 = new DNASequence(before + insertRegionProximal);
                    String ampliconName = id + "-" + jd.getFeatureLabel(END_TYPE.promimal);
                    seq1.setAccession(new AccessionID(ampliconName));
                    seq1.setDescription(idFull);
                    TextFeature tf = new TextFeature<>(jd.getFeatureLabel(END_TYPE.promimal), "Vector", jd.junctionName, jd.insertRegionProximal.toString());
                    tf.setLocation(new SequenceLocation<>( before.length() + 1, seq1.getBioEnd(), seq1, Strand.POSITIVE));
                    seq1.addFeature(tf);
                    runPrimer3(ampliconName, jd.getFeatureLabel(END_TYPE.promimal), seq1, before.length(), outputTsv.getParentFile(), primerWriter);

                    seq2 = new DNASequence(insertRegionDistal + after);
                    String ampliconName2 = id + "-" + jd.getFeatureLabel(END_TYPE.distal);
                    seq2.setAccession(new AccessionID(ampliconName2));
                    seq2.setDescription(idFull);
                    TextFeature tf2 = new TextFeature<>(jd.getFeatureLabel(END_TYPE.distal), "Vector", jd.junctionName, jd.insertRegionDistal.toString());
                    tf2.setLocation(new SequenceLocation<>(1, insertRegionDistal.length(), seq2, Strand.POSITIVE));
                    seq2.addFeature(tf2);
                    runPrimer3(ampliconName2, jd.getFeatureLabel(END_TYPE.distal), seq2, insertRegionDistal.length(), outputTsv.getParentFile(), primerWriter);
                }
            }
            catch (CompoundNotFoundException e) {
                throw new GATKException(e.getMessage(), e);
            }

            return Pair.of(seq1, seq2);
        }

        private void runPrimer3(String siteName, String junctionName, DNASequence sequence, int junctionSite, File outDir, CSVWriter primerWriter) {
            logger.info("Running primer3 for: " + siteName + ", " + junctionName);

            String seqString = sequence.getSequenceAsString();

            List<String> args = new ArrayList<>();
            args.add(primer3ExePath);

            File output = new File(outDir, "p3.output.txt");
            args.add("--output=" + output.getPath());

            args.add("--default_version=2");
            args.add("--strict_tags");

            File errorFile = new File(outDir, "p3.error.txt");
            args.add("--error=" + errorFile.getPath());

            int start = junctionSite - 50;
            File input = new File(outDir, "p3.input.txt");
            try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(input)) {
                writer.write(
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
            }
            catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            args.add(input.getPath());

            final ProcessSettings prs = new ProcessSettings(args.toArray(new String[args.size()]));
            prs.getStderrSettings().printStandard(true);
            prs.getStdoutSettings().printStandard(true);
            ProcessController.getThreadLocal().exec(prs);

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
            }
            catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }

            input.delete();
            output.delete();
            errorFile.delete();

            if (!outputMap.containsKey("PRIMER_PAIR_NUM_RETURNED")) {
                logger.info("primer3 returned no pairs");
            }

            int numPairs = Integer.parseInt(outputMap.get("PRIMER_PAIR_NUM_RETURNED"));
            for (int i=0;i<numPairs;i++) {
                String primer1Seq = outputMap.get("PRIMER_LEFT_" + i + "_SEQUENCE");
                int primer1Loc = Integer.parseInt(outputMap.get("PRIMER_LEFT_" + i).split(",")[0]) + 1;

                String primer2Seq = outputMap.get("PRIMER_RIGHT_" + i + "_SEQUENCE");
                int primer2Loc = Integer.parseInt(outputMap.get("PRIMER_RIGHT_" + i).split(",")[0]) + 1;

                String pairName = "Pair" + (i+1);
                if (primerWriter != null) {
                    primerWriter.writeNext(new String[]{siteName, junctionName, pairName, primer1Seq, String.valueOf(primer1Loc), primer2Seq, String.valueOf(primer2Loc)});
                }

                TextFeature tf = new TextFeature<>(pairName + "-F", "Primer3", pairName + "-F", pairName + "-F");
                tf.setLocation(new SequenceLocation<>(primer1Loc, primer1Loc + primer1Seq.length(), sequence, Strand.POSITIVE));
                sequence.addFeature(tf);

                TextFeature tf2 = new TextFeature<>(pairName + "-R", "Primer3", pairName + "-R", pairName + "-R");
                tf2.setLocation(new SequenceLocation<>(primer2Loc, primer2Loc + primer2Seq.length(), sequence, Strand.NEGATIVE));
                sequence.addFeature(tf2);
            }
        }
    }
}


