package com.github.discvrseq.walkers;

import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.io.FileUtils;
import org.apache.lucene.document.Document;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class VcfToLuceneIndexerIntegrationTest extends BaseIntegrationTest {

    public static String VCF4headerStrings =
                    "##fileformat=VCFv4.0\n"+
                    "##filedate=2010-06-21\n"+
                    "##reference=NCBI36\n"+
                    "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                    "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                    "##INFO=<ID=AF, Number=A, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                    "##INFO=<ID=CA, Number=1, Type=String, Description=\"Pilot 1 callability mask\">\n"+
                    "##INFO=<ID=HP, Number=1, Type=Integer, Description=\"Reference homopolymer tract length\">\n"+
                    "##INFO=<ID=NS, Number=1, Type=Integer, Description=\"Number of samples with data\">\n"+
                    "##INFO=<ID=DB, Number=0, Type=Flag, Description=\"dbSNP membership build 129 - type match and indel sequence length match within 25 bp\">\n"+
                    "##INFO=<ID=NR, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on reverse strand\">\n"+
                    "##INFO=<ID=NF, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on forward strand\">\n"+
                    "##INFO=<ID=F1, Number=1, Type=String, Description=\"Test field 1\">\n"+
                    "##INFO=<ID=F2, Number=1, Type=String, Description=\"Test field 2\">\n"+
                    "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                    "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                    "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                    "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    @Test
    public void doBasicTest() {
        Document doc = new Document();
        VCFHeader header = createHeader(VCF4headerStrings);

        final String field1Value = "Test1";
        final String field2Value = "Test2";

        Map<String, List<Object>> testFields = new HashMap<>();
        testFields.put("F1", Collections.singletonList(field1Value));
        testFields.put("F2", Collections.singletonList(field2Value));

        VcfToLuceneIndexer.addFieldsToDocument(doc, header, testFields);

        Assert.assertEquals(doc.get("F1"), field1Value);
        Assert.assertEquals(doc.get("F2"), field2Value);
    }

    private VCFHeader createHeader(String headerStr) {
        VCFCodec codec = new VCFCodec();
        return (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new ByteArrayInputStream(headerStr.getBytes(StandardCharsets.UTF_8)))));
    }

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

        args.addRaw("-AN");
        args.addRaw("SampleList");

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
    }
}
