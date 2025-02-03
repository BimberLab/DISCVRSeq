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
        args.addFlag("discardNonVariantSites");

        File sampleMappingFile = new File(outDir, "sampleMappingFile.txt");
        File outputVcf1 = new File(outDir, "outputVcf1.vcf");
        File outputVcf2 = new File(outDir, "outputVcf2.vcf");
        File outputVcf3 = new File(outDir, "outputVcf3.vcf");

        try (ICSVWriter writer = CsvUtils.getTsvWriter(sampleMappingFile)) {
            writer.writeNext(new String[]{outputVcf1.getPath(), "Sample1"});
            writer.writeNext(new String[]{outputVcf2.getPath(), "Sample1"});
            writer.writeNext(new String[]{outputVcf2.getPath(), "Sample3"});
            writer.writeNext(new String[]{outputVcf3.getPath(), "Sample3"});
        }

        args.add("sample-mapping-file", normalizePath(sampleMappingFile));

        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(outputVcf1);
        String expectedMD5 = Utils.calculateFileMD5(getTestFile(outputVcf1.getName()));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(outputVcf2);
        expectedMD5 = Utils.calculateFileMD5(getTestFile(outputVcf2.getName()));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(outputVcf3);
        expectedMD5 = Utils.calculateFileMD5(getTestFile(outputVcf3.getName()));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    private ArgumentsBuilder getBaseArgs(File outDir) {
        return getBaseArgs(outDir, 1);
    }

    @Test
    public void testInfoSubsetting() throws Exception {
        final File outDir = IOUtils.createTempDir("splitVcfBySamples.");
        ArgumentsBuilder args = getBaseArgs(outDir, 1, "mergeVcfWithAlts.vcf");
        args.addRaw("--recalculate-ac");
        args.addRaw("--keep-original-ac");
        args.addRaw("--original-ac-suffix");
        args.addRaw(".new");

        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcfWithAlts.1of2.vcf"));
        String expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcfWithAlts.1of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcfWithAlts.2of2.vcf"));
        expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcfWithAlts.2of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test
    public void testInfoSubsettingRemoveAlts() throws Exception {
        final File outDir = IOUtils.createTempDir("splitVcfBySamples.");
        ArgumentsBuilder args = getBaseArgs(outDir, 1, "mergeVcfWithAlts.vcf");
        args.addRaw("--recalculate-ac");
        args.addRaw("--keep-original-ac");
        args.addRaw("--original-ac-suffix");
        args.addRaw(".new");
        args.addRaw("--remove-unused-alternates");
        args.addFlag("discardNonVariantSites");

        runCommandLine(args);

        String actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcfWithAlts.1of2.vcf"));
        String expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcfWithAltsRemoveAlts.1of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);

        actualMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcfWithAlts.2of2.vcf"));
        expectedMD5 = Utils.calculateFileMD5(getTestFile("mergeVcfWithAltsRemoveAlts.2of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    private ArgumentsBuilder getBaseArgs(File outDir, @Nullable Integer samplesPerVcf) {
        return getBaseArgs(outDir, samplesPerVcf, "mergeVcf3.vcf");
    }

    private ArgumentsBuilder getBaseArgs(File outDir, @Nullable Integer samplesPerVcf, String inputVcf) {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, inputVcf);
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
