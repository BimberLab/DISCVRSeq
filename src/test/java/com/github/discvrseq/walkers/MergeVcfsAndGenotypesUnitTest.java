package com.github.discvrseq.walkers;

import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.StringBufferInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Test out pieces of the MergeVcfsAndGenotypes code. This is adapted from GATK3 code.
 *
 */
public class MergeVcfsAndGenotypesUnitTest {

    public static int VCF4headerStringCount = 16;

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
                    "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                    "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                    "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                    "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // this header is a small subset of the header in VCFHeaderUnitTest: VCF4headerStrings
    public static String VCF4headerStringsSmallSubset =
                "##fileformat=VCFv4.0\n" +
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    public static String VCF4headerStringsWithSamplesName =
                "##fileformat=VCFv4.0\n" +
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891\n";

    public static String VCF4headerStringsWithUniqueSamplesName =
                "##fileformat=VCFv4.0\n" +
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12892\n";


    // altered info field
    public static String VCF4headerStringsBrokenInfo =
                "##fileformat=VCFv4.0\n"+
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=1, Type=String, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+ // string to integer
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // altered format field
    public static String VCF4headerStringsBrokenFormat =
                "##fileformat=VCFv4.0\n"+
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=6, Type=String, Description=\"Genotype\">\n"+ // changed 1 to 6 here
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    private VCFHeader createHeader(String headerStr) {
        VCFCodec codec = new VCFCodec();
        VCFHeader head = null;
        head = (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(headerStr))));
        return head;
    }

    @Test
    public void testHeadersWithSamplesNamesDuplicationThatIsNotAllowed() {
        VCFHeader one = createHeader(VCF4headerStringsWithSamplesName);
        VCFHeader two = createHeader(VCF4headerStringsWithSamplesName);
        Map<String, VCFHeader> headers = new HashMap<String, VCFHeader>();
        headers.put("VCF4headerStringsWithSamplesName",one);
        headers.put("VCF4headerStringsWithSamplesName2",two);
        Assert.assertFalse(MergeVcfsAndGenotypes.areSamplesUnique(headers.values()));
    }

    @Test
    public void testHeadersWithoutSamplesNamesDuplication() {
        VCFHeader one = createHeader(VCF4headerStringsWithSamplesName);
        VCFHeader two = createHeader(VCF4headerStringsWithUniqueSamplesName);
        Map<String, VCFHeader> headers = new HashMap<String, VCFHeader>();
        headers.put("VCF4headerStringsWithSamplesName",one);
        headers.put("VCF4headerStringsWithSamplesName2",two);
        Assert.assertTrue(MergeVcfsAndGenotypes.areSamplesUnique(headers.values()));
    }

    @Test
    public void testHeadersWhereOneIsAStrictSubsetOfTheOther() {
        VCFHeader one = createHeader(VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsSmallSubset);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, false);
        Assert.assertEquals(lines.size(), VCF4headerStringCount);
    }

    @Test(expectedExceptions=IllegalStateException.class)
    public void testHeadersInfoDifferentValues() {
        VCFHeader one = createHeader(VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsBrokenInfo);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, false);
        Assert.assertEquals(lines.size(), VCF4headerStringCount);
    }

    @Test
    public void testHeadersFormatDifferentValues() {
        VCFHeader one = createHeader(VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsBrokenFormat);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, false);
        Assert.assertEquals(lines.size(), VCF4headerStringCount);
    }
}
