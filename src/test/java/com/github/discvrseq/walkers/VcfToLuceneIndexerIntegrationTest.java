package com.github.discvrseq.walkers;

import org.apache.commons.io.FileUtils;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Term;
import org.apache.lucene.queryparser.classic.MultiFieldQueryParser;
import org.apache.lucene.queryparser.classic.ParseException;
import org.apache.lucene.queryparser.flexible.standard.StandardQueryParser;
import org.apache.lucene.queryparser.flexible.standard.config.PointsConfig;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

public class VcfToLuceneIndexerIntegrationTest extends BaseIntegrationTest {

    private ArgumentsBuilder getBaseArgs(File luceneOutDir)
    {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-O");
        args.addRaw(normalizePath(luceneOutDir));

        args.addRaw("-IF");
        args.addRaw("PURPOSE");

        args.addRaw("-IF");
        args.addRaw("MAC");

        args.addRaw("-IF");
        args.addRaw("RN");

        args.addRaw(" --tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }

    @Test
    public void doBasicTest2() throws Exception {
        File luceneOutDir = new File(getTmpDir(), "luceneOutDir");
        if (luceneOutDir.exists())
        {
            FileUtils.deleteDirectory(luceneOutDir);
        }

        ArgumentsBuilder args = getBaseArgs(luceneOutDir);
        runCommandLine(args);

        File[] outputs = luceneOutDir.listFiles();
        Assert.assertEquals(5, outputs.length);

        validateLuceneIndex(luceneOutDir);
    }

    @Test
    public void doBasicTestWithMultithreading() throws Exception {
        File luceneOutDir = new File(getTmpDir(), "luceneOutDir2");
        if (luceneOutDir.exists())
        {
            FileUtils.deleteDirectory(luceneOutDir);
        }

        ArgumentsBuilder args = getBaseArgs(luceneOutDir);
        args.addRaw("--threads");
        args.addRaw("2");
        runCommandLine(args);

        File[] outputs = luceneOutDir.listFiles();
        Assert.assertEquals(5, outputs.length);

        validateLuceneIndex(luceneOutDir);
    }

    private void validateLuceneIndex(File indexPath) throws IOException, ParseException {
        try (Directory indexDirectory = FSDirectory.open(indexPath.toPath());
             IndexReader indexReader = DirectoryReader.open(indexDirectory)
        ) {
            IndexSearcher indexSearcher  = new IndexSearcher(indexReader);

            MultiFieldQueryParser queryParser = new MultiFieldQueryParser(new String[]{"contig", "start", "PURPOSE", "genomicPosition", "Samples"}, new StandardAnalyzer());

            TopDocs topDocs = indexSearcher.search(queryParser.parse("contig:=1"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 37L);

            topDocs = indexSearcher.search(new TermQuery(new Term("contig", "1")), 10);
            Assert.assertEquals(topDocs.totalHits.value, 37L);

            topDocs = indexSearcher.search(queryParser.parse("PURPOSE:=diff_pos_same_ref_same_alt"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            topDocs = indexSearcher.search(IntPoint.newRangeQuery("start", 0, 3000), 10);
            Assert.assertEquals(topDocs.totalHits.value, 10L);

            topDocs = indexSearcher.search(IntPoint.newRangeQuery("genomicPosition", 0, 3000), 10);
            Assert.assertEquals(topDocs.totalHits.value, 10L);
        }
    }

    private ArgumentsBuilder getExtendedArgs(File luceneOutDir)
    {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "indexInput.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-O");
        args.addRaw(normalizePath(luceneOutDir));

        args.addRaw("-IF");
        args.addRaw("AF");

        args.addRaw("-IF");
        args.addRaw("HaplotypeScore");

        args.addRaw("-IF");
        args.addRaw("REFFIELD");

        args.addRaw("-IF");
        args.addRaw("UB");

        args.addRaw("-IF");
        args.addRaw("FLAG");

        args.addRaw(" --tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }

    @Test
    public void doExtendedTest() throws Exception {
        File luceneOutDir = new File(getTmpDir(), "luceneOutDir");
        if (luceneOutDir.exists())
        {
            FileUtils.deleteDirectory(luceneOutDir);
        }

        ArgumentsBuilder args = getExtendedArgs(luceneOutDir);
        runCommandLine(args);

        File[] outputs = luceneOutDir.listFiles();
        Assert.assertEquals(5, outputs.length);

        try (Directory indexDirectory = FSDirectory.open(luceneOutDir.toPath());
             IndexReader indexReader = DirectoryReader.open(indexDirectory)
        ) {
            IndexSearcher indexSearcher  = new IndexSearcher(indexReader);

            MultiFieldQueryParser queryParser = new MultiFieldQueryParser(new String[]{"contig", "start", "PURPOSE", "genomicPosition", "Samples", "ref", "REFFIELD", "variableSamples"}, new StandardAnalyzer());

            StandardQueryParser numericQueryParser = new StandardQueryParser();
            numericQueryParser.setAnalyzer(new StandardAnalyzer());

            PointsConfig pointsConfig = new PointsConfig(new DecimalFormat(), Double.class);
            Map<String, PointsConfig> pointsConfigMap = new HashMap<>();
            pointsConfigMap.put("HaplotypeScore", pointsConfig);
            pointsConfigMap.put("UB", pointsConfig);
            numericQueryParser.setPointsConfigMap(pointsConfigMap);

            // Documents where contig == 1.
            TopDocs topDocs = indexSearcher.search(queryParser.parse("contig:=1"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 6L);

            // Documents where contig == 1.
            topDocs = indexSearcher.search(new TermQuery(new Term("contig", "1")), 10);
            Assert.assertEquals(topDocs.totalHits.value, 6L);

            // Documents where REFFIELD == G.
            topDocs = indexSearcher.search(queryParser.parse("REFFIELD:G"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where start >= 0, <= 65, with the rangeQuery function.
            topDocs = indexSearcher.search(IntPoint.newRangeQuery("start", 0, 65), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where HaplotypeScore == 12, with query syntax.
            topDocs = indexSearcher.search(numericQueryParser.parse("HaplotypeScore:[12.0 TO 12.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            // Documents where HaplotypeScore <= 10.0, with the query syntax.
            topDocs = indexSearcher.search(numericQueryParser.parse("HaplotypeScore:[* TO 10.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            // Documents where ref == T.
            topDocs = indexSearcher.search(queryParser.parse("ref:T"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 2L);

            // Documents where UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where UB > 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0 TO *]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where UB > 1.0 OR UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0 TO *] OR UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where UB > 1.0 AND UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0 TO *] AND UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where UB < 0.99 (should be none).
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[* TO 0.99]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 0L);

            // Documents where UB < 0.99 AND UB == 1.0 (should be none).
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[* TO 0.99] AND UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 0L);

            // Documents where UB < 0.99 AND UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[* TO 0.99] OR UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where Sample2 is variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:Sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 2L);

            // Documents where Sample2 is not variable. (All documents, minus the Sample2's.)
            topDocs = indexSearcher.search(queryParser.parse("*:* -variableSamples:Sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 4L);

            // Documents where Sample1 and Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample1 AND Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 2L);

            // Documents where Sample2 and Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample2 AND Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where Sample2 or Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample2 OR Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            // Documents where Sample1 and Sample3 are variable, and where Sample2 is NOT variable.
            // + requires that variableSamples exists, otherwise we return documents with no variableSamples field as well.
            topDocs = indexSearcher.search(queryParser.parse("*:* +variableSamples:(Sample1 AND Sample3) -variableSamples:Sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where Sample1 and Sample3 are variable, and where Sample2 is NOT variable. Tests for case insensitivity.
            topDocs = indexSearcher.search(queryParser.parse("*:* +variableSamples:(sample1 AND sample3) -variableSamples:sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Documents where Sample1, Sample2, OR Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample1 OR Sample2 OR Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            // Documents with none of (Sample1, Sample2, Sample3).
            // Said another way, at least one of Sample1, Sample2, or Sample3 is variable
            topDocs = indexSearcher.search(queryParser.parse("*:* -variableSamples:(Sample1 OR Sample2 OR Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            // Documents that don't contain one of (Sample1, Sample2, Sample3)
            topDocs = indexSearcher.search(queryParser.parse("*:* -variableSamples:Sample1 OR -variableSamples:Sample2 OR -variableSamples:Sample3"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);
        }
    }
}
