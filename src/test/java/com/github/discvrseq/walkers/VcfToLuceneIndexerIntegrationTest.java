package com.github.discvrseq.walkers;

import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.io.FileUtils;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Term;
import org.apache.lucene.queryparser.classic.MultiFieldQueryParser;
import org.apache.lucene.queryparser.classic.ParseException;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
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
        args.addRaw("REF");

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
}
