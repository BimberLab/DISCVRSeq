package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

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

    private ArgumentsBuilder getBaseArgs(File outDir) {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "mergeVcf3.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.add("R", normalizePath(getHg19Micro()));

        args.add("samplesPerVcf", 1);

        args.addRaw("-O");
        args.addRaw(normalizePath(outDir));

        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }
}
