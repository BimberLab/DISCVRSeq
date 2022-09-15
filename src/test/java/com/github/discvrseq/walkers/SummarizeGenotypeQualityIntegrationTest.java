package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class SummarizeGenotypeQualityIntegrationTest extends BaseIntegrationTest {
    @Test
    public void basicTest() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File input = new File(testBaseDir, "mergeVcf3.vcf");
        args.add("V", normalizePath(input));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(normalizePath(getTestFile("/basicTest.txt"))));

        spec.executeTest("basicTest", this);
    }
}
