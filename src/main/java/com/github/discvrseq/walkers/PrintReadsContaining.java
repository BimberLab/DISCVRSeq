package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;

import javax.annotation.Nullable;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
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

    @Argument(fullName = "output", doc = "The output file for the first FASTQ file", optional = false)
    private File outputFile1 = null;

    @Argument(fullName = "output2", doc = "The output file for the second FASTQ file.  Required if -FQ2 is used.", optional = true)
    private File outputFile2 = null;

    @Argument(fullName = "matchAllExpressions", shortName = "ma", doc = "If provided, a read pair must match all sequences/expressions to be included.  Default: false", optional = true)
    private boolean matchAllExpressions = false;

    @Argument(fullName = "expressions", shortName = "e", doc = "A list of sequences or regular expressions to test.  If either the forward or reverse read matches, it will be included.", optional = true)
    private List<String> expressions = new ArrayList<>();

    @Argument(fullName = "read1Expressions", shortName = "e1", doc = "A list of sequences or regular expressions to test against the forward read. If passing, the pair will be included.", optional = true)
    private List<String> read1Expressions = new ArrayList<>();

    @Argument(fullName = "read2Expressions", shortName = "e2", doc = "A list of sequences or regular expressions to test against the reverse read. If passing, the pair will be included.", optional = true)
    private List<String> read2Expressions = new ArrayList<>();

    private List<Pattern> eitherReadPatterns = new ArrayList<>();
    private List<Pattern> read1Patterns = new ArrayList<>();
    private List<Pattern> read2Patterns = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsWritable(outputFile1);

        if (FASTQ2 != null) {
            IOUtil.assertFileIsReadable(FASTQ2);
            IOUtil.assertFileIsWritable(outputFile2);
        }

        if (!read2Expressions.isEmpty() && FASTQ2 == null){
            throw new UserException.BadInput("Specified --read2Expressions, but --fastq2 was not provided");
        }

        //initialize expressions:
        expressions.forEach(e -> {
            eitherReadPatterns.add(Pattern.compile(e));
        });

        read1Expressions.forEach(e -> {
            read1Patterns.add(Pattern.compile(e));
        });

        read2Expressions.forEach(e -> {
            read2Patterns.add(Pattern.compile(e));
        });
    }

    @Override
    public void traverse() {
        FastqWriterFactory fact = new FastqWriterFactory();
        fact.setUseAsyncIo(true);

        long written = 0L;
        try (FastqReader reader1 = fileToFastqReader(FASTQ); FastqReader reader2 = FASTQ2 == null ? null : fileToFastqReader(FASTQ2);FastqWriter writer1 = fact.newWriter(outputFile1); FastqWriter writer2 = FASTQ2 == null ? null : fact.newWriter(outputFile2)) {
            while(reader1.hasNext())
            {
                FastqRecord fq1 = reader1.next();
                FastqRecord fq2 = reader2 == null ? null : reader2.next();

                if (acceptReadPair(fq1, fq2)) {
                    writer1.write(fq1);
                    if (writer2 != null){
                        writer2.write(fq2);
                    }
                    written++;
                }
            }
        }

        logger.info("total reads written: " + written);
    }

    private FastqReader fileToFastqReader(final File file) {
        return new FastqReader(file, true);
    }

    public boolean acceptReadPair(FastqRecord read1, @Nullable FastqRecord read2) {
        if (!eitherReadPatterns.isEmpty()) {
            int matchesPair = inspect(eitherReadPatterns, read1, read2);
            if (!isPassing(matchesPair, eitherReadPatterns)) {
                return false;
            }
        }

        if (!read1Patterns.isEmpty()) {
            int matches1 = inspect(read1Patterns, read1);
            if (!isPassing(matches1, read1Patterns)) {
                return false;
            }
        }

        if (!read2Patterns.isEmpty()) {
            if (read2 == null) {
                throw new UserException.BadInput("Specified read2 expressions, but read2 not found");
            }
            else {
                int matches2 = inspect(read2Patterns, read2);
                if (!isPassing(matches2, read2Patterns)) {
                    return false;
                }
            }
        }

        return true;
    }

    private boolean isPassing(int matches, List<Pattern> expressions) {
        return matchAllExpressions ? matches == expressions.size() : matches > 0;
    }

    private int inspect(List<Pattern> exprs, FastqRecord... reads) {
        int matching = 0;
        for (Pattern expr : exprs) {
            for (FastqRecord read : reads) {
                Matcher m = expr.matcher(read.getReadString());
                if (m.find()) {
                    matching++;
                    break;
                }
            }
        }

        return matching;
    }
}
