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
        ArgumentsBuilder args = getBaseArgs();

        final File outDir = IOUtils.createTempDir("splitVcfBySamples.");

        args.addRaw("-O");
        args.addRaw(normalizePath(outDir));

        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        runCommandLine(args);

        String expectedMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.1of2.vcf"));
        String actualMD5 = Utils.calculateFileMD5(getTestFile("mergeVcf3.1of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);

        expectedMD5 = Utils.calculateFileMD5(new File(outDir, "mergeVcf3.2of2.vcf"));
        actualMD5 = Utils.calculateFileMD5(getTestFile("mergeVcf3.2of2.vcf"));
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "mergeVcf3.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.add("R", normalizePath(getHg19Micro()));

        args.add("samplesPerVcf", 1);

        return args;
    }
}
