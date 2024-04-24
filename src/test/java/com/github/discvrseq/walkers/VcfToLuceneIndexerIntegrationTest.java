package com.github.discvrseq.walkers;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.commons.io.FileUtils;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Term;
import org.apache.lucene.queryparser.classic.MultiFieldQueryParser;
import org.apache.lucene.queryparser.classic.ParseException;
import org.apache.lucene.queryparser.flexible.standard.StandardQueryParser;
import org.apache.lucene.queryparser.flexible.standard.config.PointsConfig;
import org.apache.lucene.search.*;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.net.URI;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
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

        args.addRaw("--index-stats");
        args.addRaw(new File(luceneOutDir.getPath() + ".stats.txt"));

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

        // Validate stats:
        File expected = getTestFile("luceneOutDir.stats.txt");
        File actual = new File(luceneOutDir.getPath() + ".stats.txt");
        IntegrationTestSpec.assertEqualTextFiles(actual, expected);

        // Validate index
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

            ScoreDoc[] scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                String fieldValue = document.get("contig");

                Assert.assertEquals("1", fieldValue);
            }

            // Documents where REFFIELD == G.
            topDocs = indexSearcher.search(queryParser.parse("REFFIELD:G"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                String fieldValue = document.get("REFFIELD");

                Assert.assertEquals("G", fieldValue);
            }

            // Documents where start >= 0, <= 65, with the rangeQuery function.
            topDocs = indexSearcher.search(IntPoint.newRangeQuery("start", 0, 65), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                int fieldValue = Integer.parseInt(document.get("start"));

                Assert.assertTrue(fieldValue >= 0);
                Assert.assertTrue(fieldValue <= 65);
            }

            // Documents where HaplotypeScore == 12, with query syntax.
            topDocs = indexSearcher.search(numericQueryParser.parse("HaplotypeScore:[12.0 TO 12.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                float fieldValue = Float.parseFloat(document.get("HaplotypeScore"));

                Assert.assertEquals(12.0, fieldValue);
            }

            // Documents where HaplotypeScore <= 10.0, with the query syntax.
            topDocs = indexSearcher.search(numericQueryParser.parse("HaplotypeScore:[* TO 10.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                Float fieldValue = Float.parseFloat(document.get("HaplotypeScore"));

                Assert.assertTrue(fieldValue <= 10.0);
            }

            // Documents where ref == T.
            topDocs = indexSearcher.search(queryParser.parse("ref:T"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 2L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                String fieldValue = document.get("ref");

                Assert.assertEquals("T", fieldValue);
            }

            // Documents where UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                float fieldValue = Float.parseFloat(document.get("UB"));

                Assert.assertEquals(1.0, fieldValue);
            }

            // Documents where UB >= 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0 TO *]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                float fieldValue = Float.parseFloat(document.get("UB"));

                Assert.assertTrue(fieldValue >= 1.0);
            }

            // Documents where UB > 1.0 OR UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[1.0000001 TO *] OR UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                float fieldValue = Float.parseFloat(document.get("UB"));

                Assert.assertTrue(fieldValue > 1.0  || fieldValue == 1.0);
            }

            // Documents where UB < 0.99 (should be none).
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[* TO 0.99]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 0L);

            // Documents where UB < 0.99 AND UB == 1.0 (should be none).
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[* TO 0.99] AND UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 0L);

            // Documents where UB < 0.99 OR UB == 1.0.
            topDocs = indexSearcher.search(numericQueryParser.parse("UB:[* TO 0.99] OR UB:[1.0 TO 1.0]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                float fieldValue = Float.parseFloat(document.get("UB"));

                Assert.assertTrue(fieldValue < 0.99 || fieldValue == 1.0);
            }

            // Documents where Sample2 is variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:Sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 2L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample2"));
            }

            // Documents where Sample2 is not variable. (All documents, minus the Sample2's.)
            topDocs = indexSearcher.search(queryParser.parse("*:* -variableSamples:Sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 4L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(!variableSamples.contains("Sample2"));
            }

            // Documents where Sample1 and Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample1 AND Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 2L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample1") && variableSamples.contains("Sample3"));
            }

            // Documents where Sample2 and Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample2 AND Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample2") && variableSamples.contains("Sample3"));
            }

            // Documents where Sample2 or Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample2 OR Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample2") || variableSamples.contains("Sample3"));
            }

            // Documents where Sample1 and Sample3 are variable, and where Sample2 is NOT variable.
            // + requires that variableSamples exists, otherwise we return documents with no variableSamples field as well.
            topDocs = indexSearcher.search(queryParser.parse("*:* +variableSamples:(Sample1 AND Sample3) -variableSamples:Sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample1") && variableSamples.contains("Sample3") && !variableSamples.contains("Sample2"));
            }

            // Documents where Sample1 and Sample3 are variable, and where Sample2 is NOT variable. Tests for case insensitivity.
            topDocs = indexSearcher.search(queryParser.parse("*:* +variableSamples:(sample1 AND sample3) -variableSamples:sample2"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample1") && variableSamples.contains("Sample3") && !variableSamples.contains("Sample2"));
            }

            // Documents where Sample1, Sample2, OR Sample3 are variable.
            topDocs = indexSearcher.search(queryParser.parse("variableSamples:(Sample1 OR Sample2 OR Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(variableSamples.contains("Sample1") || variableSamples.contains("Sample3") || variableSamples.contains("Sample2"));
            }

            // Documents with none of (Sample1, Sample2, Sample3).
            topDocs = indexSearcher.search(queryParser.parse("*:* -variableSamples:(Sample1 OR Sample2 OR Sample3)"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(!variableSamples.contains("Sample1") && !variableSamples.contains("Sample3") && !variableSamples.contains("Sample2"));
            }

            // Documents that don't contain one of (Sample1, Sample2, Sample3)
            topDocs = indexSearcher.search(queryParser.parse("*:* -variableSamples:Sample1 OR -variableSamples:Sample2 OR -variableSamples:Sample3"), 10);
            Assert.assertEquals(topDocs.totalHits.value, 3L);

            scoreDocs = topDocs.scoreDocs;
            for (ScoreDoc scoreDoc : scoreDocs) {
                int docId = scoreDoc.doc;
                Document document = indexSearcher.storedFields().document(docId);
                List<String> variableSamples = Arrays.asList(document.getValues("variableSamples"));

                Assert.assertTrue(!variableSamples.contains("Sample1") || !variableSamples.contains("Sample3") || !variableSamples.contains("Sample2"));
            }

            // Documents where HaplotypeScore == 0.12, with query syntax.
            topDocs = indexSearcher.search(numericQueryParser.parse("HaplotypeScore:[0.12 TO 0.12]", ""), 10);
            Assert.assertEquals(topDocs.totalHits.value, 1L);

            // Top 50 hits are sorted by genomicPosition
            topDocs = indexSearcher.search(new MatchAllDocsQuery(), 6, new Sort(new SortField("genomicPosition", SortField.Type.INT)));
            Assert.assertEquals(6, topDocs.scoreDocs.length);

            int lastGenomicPosition = -1;
            for (ScoreDoc scoreDoc : topDocs.scoreDocs) {
                Document document = indexSearcher.storedFields().document(scoreDoc.doc);
                int currentGenomicPosition = Integer.parseInt(document.get("genomicPosition"));
                if (lastGenomicPosition != -1) {
                    Assert.assertTrue(lastGenomicPosition <= currentGenomicPosition);
                }
                lastGenomicPosition = currentGenomicPosition;
            }

            // Results are sorted by REFFIELD
            topDocs = indexSearcher.search(new MatchAllDocsQuery(), 6, new Sort(new SortField("REFFIELD", SortField.Type.STRING)));
            Assert.assertEquals(6, topDocs.scoreDocs.length);

            String lastRefField = null;
            for (ScoreDoc scoreDoc : topDocs.scoreDocs) {
                Document document = indexSearcher.storedFields().document(scoreDoc.doc);
                String currentRefField = document.get("REFFIELD");
                if (lastRefField != null) {
                    Assert.assertTrue(lastRefField.compareTo(currentRefField) <= 0);
                }
                lastRefField = currentRefField;
            }

            // Results are sorted by start
            topDocs = indexSearcher.search(new MatchAllDocsQuery(), 6, new Sort(new SortField("start", SortField.Type.INT)));
            Assert.assertEquals(6, topDocs.scoreDocs.length);

            int lastStart = -1;
            for (ScoreDoc scoreDoc : topDocs.scoreDocs) {
                Document document = indexSearcher.storedFields().document(scoreDoc.doc);
                int currentStart = Integer.parseInt(document.get("start"));
                Assert.assertTrue(lastStart <= currentStart);
                lastStart = currentStart;
            }

            // Results are sorted by HaplotypeScore
            topDocs = indexSearcher.search(new MatchAllDocsQuery(), 6, new Sort(new SortField("HaplotypeScore", SortField.Type.DOUBLE)));
            Assert.assertEquals(6, topDocs.scoreDocs.length);

            float lastHaplotypeScore = -1.0f;
            for (ScoreDoc scoreDoc : topDocs.scoreDocs) {
                Document document = indexSearcher.storedFields().document(scoreDoc.doc);
                float currentHaplotypeScore = Float.parseFloat(document.get("HaplotypeScore"));
                Assert.assertTrue(lastHaplotypeScore <= currentHaplotypeScore);
                lastHaplotypeScore = currentHaplotypeScore;
            }

            // Results are sorted by genomicPosition
            topDocs = indexSearcher.search(new MatchAllDocsQuery(), 6, Sort.INDEXORDER);
            Assert.assertTrue(topDocs.scoreDocs.length > 0);

            lastGenomicPosition = -1;
            for (ScoreDoc scoreDoc : topDocs.scoreDocs) {
                Document document = indexSearcher.storedFields().document(scoreDoc.doc);
                int currentGenomicPosition = Integer.parseInt(document.get("genomicPosition"));
                if (lastGenomicPosition != -1) {
                    Assert.assertTrue(lastGenomicPosition <= currentGenomicPosition);
                }
                lastGenomicPosition = currentGenomicPosition;
            }
        }
    }

    @Test
    public void validateLuceneVersions() throws Exception
    {
        // One of the main use cases of lucene is to build VCF indexes read by LabKey's code. This check tries to ensure
        // the lucene versions are in sync between the two
        String discvrVersion = null;
        try (BufferedReader reader = new BufferedReader(new FileReader(new File("build.gradle")))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.contains("lucene.version")) {
                    discvrVersion = line.split("'")[3];
                    break;
                }
            }
        }

        Assert.assertNotNull(discvrVersion, "Unable to find labkey lucene version");

        String branch;
        try (CloseableHttpClient httpClient = HttpClients.createDefault()) {
            HttpGet post = new HttpGet(URI.create("https://api.github.com/repos/bimberlab/DiscvrLabKeyModules"));
            post.setHeader("Accept", "application/vnd.github+json");
            post.setHeader("X-GitHub-Api-Version", "2022-11-28");
            post.setHeader("Content-type", "application/json");

            HttpResponse response = httpClient.execute(post);
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(response.getEntity().getContent()))) {
                JsonElement json = JsonParser.parseReader(reader);
                JsonObject obj = json.getAsJsonObject();
                branch = obj.get("default_branch").getAsString();
                branch = branch.replaceAll("discvr-", "release") + "-SNAPSHOT";
            }
        }
        catch (IOException e) {
            throw e;
        }

        Assert.assertNotNull(branch, "Missing branch");

        String labkeyLucene = null;
        try (CloseableHttpClient httpClient = HttpClients.createDefault()) {
            HttpGet post = new HttpGet(URI.create("https://raw.githubusercontent.com/labkey/server/" + branch + "/gradle.properties"));
            post.setHeader("Accept", "application/vnd.github+json");
            post.setHeader("X-GitHub-Api-Version", "2022-11-28");
            post.setHeader("Content-type", "application/json");

            HttpResponse response = httpClient.execute(post);
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(response.getEntity().getContent()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("luceneVersion")) {
                        labkeyLucene = line.split("=")[1];
                        break;
                    }
                }
            }
        }
        catch (IOException e) {
            throw e;
        }

        Assert.assertNotNull(labkeyLucene, "Unable to find labkey lucene version");
        Assert.assertEquals(discvrVersion, labkeyLucene, "DISCVR-seq and LabKey Lucene versions are not equal");
    }
}
