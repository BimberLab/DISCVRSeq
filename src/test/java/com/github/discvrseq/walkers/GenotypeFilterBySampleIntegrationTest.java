package com.github.discvrseq.walkers;

import com.github.discvrseq.Main;
import htsjdk.tribble.bed.BEDCodec;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class GenotypeFilterBySampleIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        String fn = "GenotypeFilterBySample.vcf";

        args.add("-O");
        File outFile = new File(normalizePath(getSafeNonExistentFile(fn)));
        args.add(outFile);

        runCommandLine(args.getArgsArray());

        File expected = getTestFile(fn);
        IntegrationTestSpec.assertEqualTextFiles(outFile, expected);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File testBaseDir = new File(publicTestDir + "com/github/discvrseq/TestData");

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
