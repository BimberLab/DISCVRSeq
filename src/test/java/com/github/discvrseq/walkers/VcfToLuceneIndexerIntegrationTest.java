package com.github.discvrseq.walkers;

import org.apache.lucene.document.Document;
import org.apache.lucene.index.IndexableField;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import org.testng.Assert;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    public void doBasicTest() throws Exception {
        Document doc = new Document();
        VCFHeader header = createHeader(VCF4headerStrings);

        Map<String, List<Object>> testFields = new HashMap<String, List<Object>>();
        testFields.put("F1", Arrays.asList("Test1"));
        testFields.put("F2", Arrays.asList("Test2"));

        String field1Value = "Test1";
        String field2Value = "Test2";

        VcfToLuceneIndexer.addFieldsToDocument(doc, header, testFields);

        Assert.assertEquals(doc.get("F1"), field1Value);
        Assert.assertEquals(doc.get("F2"), field2Value);
    }

    private VCFHeader createHeader(String headerStr) {
        VCFCodec codec = new VCFCodec();
        VCFHeader head = null;
        head = (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new ByteArrayInputStream(headerStr.getBytes(StandardCharsets.UTF_8)))));
        return head;
    }
}
