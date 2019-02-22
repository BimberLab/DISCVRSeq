package com.github.discvrseq.walkers;

import htsjdk.tribble.bed.BEDCodec;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class GenotypeFilterBySampleIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        String fn = "GenotypeFilterBySample.vcf";

        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest("testBasicOperation", this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);

        args.add(normalizePath(input));

        args.add("-bl");
        File blackListBed = new File(testBaseDir, "GenotypeFilterBySampleBlacklist.bed");
        ensureIndex(blackListBed, new BEDCodec());
        args.add(normalizePath(blackListBed));

        return args;
    }
}
