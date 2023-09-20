package com.github.discvrseq.walkers;

import com.github.discvrseq.util.CsvUtils;
import com.opencsv.CSVWriter;
import com.opencsv.ICSVWriter;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import javax.annotation.Nullable;
import java.io.File;

public class SplitVcfBySamplesIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        final File outDir = IOUtils.createTempDir("splitVcfBySamples.");
        ArgumentsBuilder args = getBaseArgs(outDir);
        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.1of2.vcf"));
        String expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcf3.1of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.2of2.vcf"));
        expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcf3.2of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test
    public void testBasicOperationDiscard() throws Exception {
        final File outDir = IOUtils.createTempDir("splitVcfBySamples.");
        ArgumentsBuilder args = getBaseArgs(outDir);
        args.addFlag("discardNonVariantSites");

        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.1of2.vcf"));
        String expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcf3.1of2Discard.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.2of2.vcf"));
        expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcf3.2of2Discard.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test
    public void testBasicOperationMinSamples() throws Exception {
        final File outDir = IOUtils.createTempDir("splitVcfBySamples.");
        ArgumentsBuilder args = getBaseArgs(outDir);
        args.add("minAllowableInFinalVcf", 2);

        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.1of1.vcf"));
        String expectedMD5 = Utils.calculateFileMD5(new File(testBaseDir, "mergeVcf3.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test
    public void testSampleCsv() throws Exception {
        final File outDir = IOUtils.createTempDir("sampleCsv.");
        ArgumentsBuilder args = getBaseArgs(outDir, 0);

        File sampleMappingFile = new File(outDir, "sampleMappingFile.txt");
        File outputVcf1 = new File(outDir, "outputVcf1.vcf");
        File outputVcf2 = new File(outDir, "outputVcf2.vcf");;

        try (ICSVWriter writer = CsvUtils.getTsvWriter(sampleMappingFile)) {
            writer.writeNext(new String[]{outputVcf1.getPath(), "Sample1"});
            writer.writeNext(new String[]{outputVcf2.getPath(), "Sample1"});
            writer.writeNext(new String[]{outputVcf2.getPath(), "Sample3"});
        }

        args.add("sampleMappingFile", normalizePath(sampleMappingFile));

        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(outputVcf1);
        String expectedMD5 = Utils.calculateFileMD5(getTestFile(outputVcf1.getName()));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(outputVcf2);
        expectedMD5 = Utils.calculateFileMD5(getTestFile(outputVcf2.getName()));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    private ArgumentsBuilder getBaseArgs(File outDir) {
        return getBaseArgs(outDir, 1);
    }

    private ArgumentsBuilder getBaseArgs(File outDir, @Nullable Integer samplesPerVcf) {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "mergeVcf3.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.add("R", normalizePath(getHg19Micro()));

        if (samplesPerVcf != null) {
            args.add("samplesPerVcf", samplesPerVcf);
        }

        args.addRaw("-O");
        args.addRaw(normalizePath(outDir));

        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }
}
