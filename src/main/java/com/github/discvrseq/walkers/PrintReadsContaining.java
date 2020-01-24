package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.text.similarity.LevenshteinDistance;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This walker will iterate FASTQ or pair of FASTQs and print any reads matching the supplied expressions.  Expressions can be simple strings or a java regular expression.
 *
 *
 * <h3>Usage examples:</h3>
 * <h4>Simple search for a sequence in either forward or reverse read:</h4>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadsContaining \
 *     -FQ fastq_R1.fastq.gz \
 *     -FQ2 fastq_R2.fastq.gz \
 *     -e 'TACG' \
 *     -O output_R1.fastq.gz \
 *     -O2 output_R2.fastq.gz
 * </pre>
 * <h4>RegEx search for reads where either forward or reverse contains only A/T</h4>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadsContaining \
 *     -FQ fastq_R1.fastq.gz \
 *     -FQ2 fastq_R2.fastq.gz \
 *     -e '^[AT]+$ \
 *     -O output_R1.fastq.gz \
 *     -O2 output_R2.fastq.gz
 * </pre>
 * <h4>Same as above, except tests forward read only</h4>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadsContaining \
 *     -FQ fastq_R1.fastq.gz \
 *     -FQ2 fastq_R2.fastq.gz \
 *     -e1 '^[AT]+$' \
 *     -O output_R1.fastq.gz \
 *     -O2 output_R2.fastq.gz
 * </pre>
 * <h4>In this example, the forward read must be exclusively A/T, and either the forward or reverse must contain AA</h4>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadsContaining \
 *     -FQ fastq_R1.fastq.gz \
 *     -FQ2 fastq_R2.fastq.gz \
 *     -e1 '^[AT]+$' \
 *     -e 'AA' \
 *     -O output_R1.fastq.gz \
 *     -O2 output_R2.fastq.gz
 * </pre>
 * <h4>In this example, one of these two conditions must be true: either read is exclusively A/T or either read contains AA</h4>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadsContaining \
 *     -FQ fastq_R1.fastq.gz \
 *     -FQ2 fastq_R2.fastq.gz \
 *     -e '^[AT]+$' \
 *     -e 'AA' \
 *     -O output_R1.fastq.gz \
 *     -O2 output_R2.fastq.gz
 * </pre>
 * <h4>Similar to the above, except both conditions must be met.</h4>
 * <pre>
 *  java -jar DISCVRseq.jar PrintReadsContaining \
 *     -FQ fastq_R1.fastq.gz \
 *     -FQ2 fastq_R2.fastq.gz \
 *     -e '^[AT]+$' \
 *     -e 'AA' \
 *     --matchAllExpressions \
 *     -O output_R1.fastq.gz \
 *     -O2 output_R2.fastq.gz
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool iterates a set of FASTQs, printing any pair that matches the supplied sequences/expressions.",
        oneLineSummary = "Print read pairs matching a list of sequences/expressions",
        programGroup = DiscvrSeqProgramGroup.class
)
public class PrintReadsContaining extends GATKTool {
    @Argument(fullName="fastq", doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data.")
    public File FASTQ;

    @Argument(fullName="fastq2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.", optional=true)
    public File FASTQ2;

    @Argument(fullName="summaryFile", doc="If provided, a TSV summary of matches will be written here.", optional=true)
    public File SUMMARY_FILE;

    @Argument(fullName = "output", doc = "The output file for the first FASTQ file", optional = false)
    private File outputFile1 = null;

    @Argument(fullName = "output2", doc = "The output file for the second FASTQ file.  Required if -FQ2 is used.", optional = true)
    private File outputFile2 = null;

    @Argument(fullName = "matchAllExpressions", shortName = "ma", doc = "If provided, a read pair must match all sequences/expressions to be included.  Default: false", optional = true)
    private boolean matchAllExpressions = false;

    @Argument(fullName = "expressions", shortName = "e", doc = "A list of sequences or regular expressions to test.  If either the forward or reverse read matches, it will be included.", optional = true)
    private List<String> expressions = new ArrayList<>();

    @Argument(fullName = "expressionNames", shortName = "en", doc = "A mechanism to allow names for each expression.  If used, the number of expression names (-en) must be equal to the number of expressions (-e).  If not provided, the expression itself will be used.", optional = true)
    private List<String> expressionNames = new ArrayList<>();

    @Argument(fullName = "read1Expressions", shortName = "e1", doc = "A list of sequences or regular expressions to test against the forward read. If passing, the pair will be included.", optional = true)
    private List<String> read1Expressions = new ArrayList<>();

    @Argument(fullName = "read1ExpressionNames", shortName = "e1n", doc = "A mechanism to allow names for each read1 expression.  If used, the number of read1 expression names (-e1n) must be equal to the number of read1 expressions (-e1).  If not provided, the expression itself will be used.", optional = true)
    private List<String> read1ExpressionNames = new ArrayList<>();

    @Argument(fullName = "read2Expressions", shortName = "e2", doc = "A list of sequences or regular expressions to test against the reverse read. If passing, the pair will be included.", optional = true)
    private List<String> read2Expressions = new ArrayList<>();

    @Argument(fullName = "read2ExpressionNames", shortName = "e2n", doc = "A mechanism to allow names for each read2 expression.  If used, the number of read2 expression names (-e2n) must be equal to the number of read2 expressions (-e2).  If not provided, the expression itself will be used.", optional = true)
    private List<String> read2ExpressionNames = new ArrayList<>();

    @Argument(fullName = "editDistance", shortName = "ed", doc = "If provided, all expressions will be treated as simple strings, and must be ATGC characters.  The tool will scan each motif against the read(s) and report a match if the edit distance is less than or equal to this threshold.", optional = true)
    private int editDistance = 0;

    private List<SeqPattern> eitherReadPatterns = new ArrayList<>();
    private List<SeqPattern> read1Patterns = new ArrayList<>();
    private List<SeqPattern> read2Patterns = new ArrayList<>();

    public class SeqPattern {
        Pattern pattern;
        String name;

        public SeqPattern(Pattern pattern, String name)
        {
            this.pattern = pattern;
            this.name = name;
        }
    }

    public class SeqMatch {
        int start0;
        int end0;
        String exprName;
        ReadType rt;

        public SeqMatch(Matcher m, String name, ReadType rt) {
            this.start0 = m.start();
            this.end0 = m.end();

            this.exprName = name;
            this.rt = rt;
        }

        //0-based indexes
        public SeqMatch(int start, int end, String name, ReadType rt) {
            this.start0 = start;
            this.end0 = end;

            this.exprName = name;
            this.rt = rt;
        }
    }

    public class SeqPairMatch {
        Collection<SeqMatch> matches;

        public SeqPairMatch(Collection<SeqMatch> matches) {
            this.matches = matches;
        }

        public Set<String> getUniqueHitNames(ReadType rt) {
            Set<String> ret = new HashSet<>();
            for (SeqMatch m : matches) {
                if (m.rt == rt || rt == ReadType.Any)
                    ret.add(m.exprName);
            }

            return ret;
        }
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsWritable(outputFile1);

        if (FASTQ2 != null) {
            IOUtil.assertFileIsReadable(FASTQ2);
            IOUtil.assertFileIsWritable(outputFile2);
        }

        if (SUMMARY_FILE != null) {
            IOUtil.assertFileIsWritable(SUMMARY_FILE);
        }

        if (!read2Expressions.isEmpty() && FASTQ2 == null){
            throw new UserException.BadInput("Specified --read2Expressions, but --fastq2 was not provided");
        }

        if (!expressionNames.isEmpty() && expressionNames.size() != expressions.size()){
            throw new UserException.BadInput("Expression names provided, but the number of expressions doesnt match the number of names (see -e and -en");
        }

        if (!read1ExpressionNames.isEmpty() && read1ExpressionNames.size() != read1Expressions.size()){
            throw new UserException.BadInput("Read1 expression names provided, but the number of read1 expressions doesnt match the number of names (see -e1 and -e1n");
        }

        if (!read2ExpressionNames.isEmpty() && read2ExpressionNames.size() != read2Expressions.size()){
            throw new UserException.BadInput("Read2 expression names provided, but the number of read2 expressions doesnt match the number of names (see -e2 and -e2n");
        }

        //initialize expressions:
        initializeExpressions(expressions, expressionNames, eitherReadPatterns);
        initializeExpressions(read1Expressions, read1ExpressionNames, read1Patterns);
        initializeExpressions(read2Expressions, read2ExpressionNames, read2Patterns);
    }

    final Pattern ntMatch = Pattern.compile("^[ATGC]+$");

    private void initializeExpressions(List<String> expressions, List<String> names, List<SeqPattern> target) {
        for (int i=0;i<expressions.size();i++) {
            String name = names.isEmpty() ? expressions.get(i) : names.get(i);
            String expr = expressions.get(i);
            //validate queries are ATGC only:
            if (editDistance > 0) {
                if (!ntMatch.matcher(expr).matches()) {
                    throw new UserException.BadInput("Expression must used only bases (ATGC) when edit distance is used: " + expr);
                }

                expr = expr.toUpperCase();
            }

            target.add(new SeqPattern(Pattern.compile(expr), name));
        };
    }

    @Override
    public void traverse() {
        FastqWriterFactory fact = new FastqWriterFactory();
        fact.setUseAsyncIo(true);

        Map<String, Long> matchCount = new HashMap<>();
        Map<String, Long> matchCountR1 = new HashMap<>();
        Map<String, Long> matchCountR2 = new HashMap<>();

        long totalReads = 0L;
        long written = 0L;
        try (FastqReader reader1 = fileToFastqReader(FASTQ); FastqReader reader2 = FASTQ2 == null ? null : fileToFastqReader(FASTQ2); FastqWriter writer1 = fact.newWriter(outputFile1); FastqWriter writer2 = FASTQ2 == null ? null : fact.newWriter(outputFile2); CSVWriter csvWriter = SUMMARY_FILE == null ? null : new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(SUMMARY_FILE), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            if (csvWriter != null) {
                csvWriter.writeNext(new String[]{"ReadName", "ReadType", "ExpressionName", "Start", "End", "TotalHitsForPair"});
            }

            while(reader1.hasNext())
            {
                FastqRecord fq1 = reader1.next();
                FastqRecord fq2 = reader2 == null ? null : reader2.next();
                totalReads++;

                SeqPairMatch matches = findMatches(fq1, fq2);
                if (matches != null) {
                    writer1.write(fq1);
                    if (writer2 != null){
                        writer2.write(fq2);
                    }
                    written++;

                    appendCounts(matchCount, matches, ReadType.Any);
                    appendCounts(matchCountR1, matches, ReadType.Forward);
                    appendCounts(matchCountR2, matches, ReadType.Reverse);

                    if (csvWriter != null) {
                        writeMatchSummary(matches, fq1, fq2, csvWriter);
                    }
                }
            }
        }
        catch (IOException e) {
            throw new GATKException("There was an error writing data", e);
        }

        logger.info("total reads inspected: " + totalReads);
        logger.info("total reads accepted: " + written);
        logger.info("the following counts were identified per expression.  note: each read pair can match multiple expressions, and these values represent the total matches, not total reads that were matched:");
        for (String name : matchCount.keySet()) {
            String perBase = " (R1: " + matchCountR1.getOrDefault(name, 0L) + " / R2: " + matchCountR2.getOrDefault(name, 0L) + ")";
            logger.info(name + ": " + matchCount.get(name) + perBase);
        }
    }

    private void writeMatchSummary(SeqPairMatch matches, FastqRecord fq1, FastqRecord fq2, CSVWriter csvWriter) {
        for (SeqMatch m : matches.matches) {
            String readType = m.rt.name();
            csvWriter.writeNext(new String[]{fq1.getReadName(), readType, m.exprName, String.valueOf(m.start0), String.valueOf(m.end0), String.valueOf(matches.matches.size())});
        }
    }

    private FastqReader fileToFastqReader(final File file) {
        return new FastqReader(file, true);
    }

    private void appendCounts(Map<String, Long> counts, SeqPairMatch match, ReadType rt) {
        for (String name : match.getUniqueHitNames(rt)){
            long val = counts.getOrDefault(name, 0L);
            val++;
            counts.put(name, val);
        }
    }

    public SeqPairMatch findMatches(FastqRecord read1, @Nullable FastqRecord read2) {
        Set<SeqMatch> matches = new HashSet<>();

        if (!eitherReadPatterns.isEmpty()) {
            Map<ReadType, FastqRecord> map = new HashMap<>();
            map.put(ReadType.Forward, read1);
            map.put(ReadType.Reverse, read2);
            List<SeqMatch> matchesPair = inspect(eitherReadPatterns, map);
            if (!isPassing(matchesPair, eitherReadPatterns)) {
                //NOTE: even if this fails, we might want to inspect the other expressions:
                if (matchAllExpressions) {
                    return null;
                }
            }
            else {
                matches.addAll(matchesPair);
            }
        }

        if (!read1Patterns.isEmpty()) {
            Map<ReadType, FastqRecord> map = new HashMap<>();
            map.put(ReadType.Forward, read1);
            List<SeqMatch> matches1 = inspect(read1Patterns, map);
            if (!isPassing(matches1, read1Patterns)) {
                if (matchAllExpressions) {
                    return null;
                }
            }
            else {
                matches.addAll(matches1);
            }
        }

        if (!read2Patterns.isEmpty()) {
            if (read2 == null) {
                throw new UserException.BadInput("Specified read2 expressions, but read2 not found");
            }
            else {
                Map<ReadType, FastqRecord> map = new HashMap<>();
                map.put(ReadType.Reverse, read2);
                List<SeqMatch> matches2 = inspect(read2Patterns, map);
                if (!isPassing(matches2, read2Patterns)) {
                    if (matchAllExpressions) {
                        return null;
                    }
                }
                else {
                    matches.addAll(matches2);
                }
            }
        }

        return matches.isEmpty() ? null : new SeqPairMatch(matches);
    }

    private boolean isPassing(List<SeqMatch> matches, List<SeqPattern> expressions) {
        return matchAllExpressions ? matches.size() == expressions.size() : !matches.isEmpty();
    }

    private List<SeqMatch> inspect(List<SeqPattern> exprs, Map<ReadType, FastqRecord> reads) {
        List<SeqMatch> matching = new ArrayList<>();
        for (SeqPattern expr : exprs) {
            for (ReadType rt : reads.keySet()) {
                if (editDistance > 0) {
                    SeqMatch match = fuzzyMatch(expr.pattern.pattern(), reads.get(rt), expr, rt);
                    if (match != null) {
                        matching.add(match);
                    }
                }
                else {
                    Matcher m = expr.pattern.matcher(reads.get(rt).getReadString());
                    if (m.find()) {
                        matching.add(new SeqMatch(m, expr.name, rt));
                        //break;  //NOTE: inspect both reads in case we have more than one hit
                    }
                }
            }
        }

        return matching;
    }

    final LevenshteinDistance levenshteinDistance = LevenshteinDistance.getDefaultInstance();

    private SeqMatch fuzzyMatch(String query, FastqRecord read, SeqPattern expr, ReadType rt) {
        String readSeq = read.getReadString().toUpperCase();
        int windows = readSeq.length() - query.length();

        int i = 0;
        while (i < windows) {
            CharSequence test = readSeq.subSequence(i, i + query.length());
            if (levenshteinDistance.apply(query, test) <= editDistance) {
                return new SeqMatch(i, query.length(), expr.name, rt);
            }

            i++;
        }

        return null;
    }

    private static enum ReadType {
        Forward(),
        Reverse(),
        Any();
    }
}
