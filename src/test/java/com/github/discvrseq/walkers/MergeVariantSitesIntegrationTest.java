package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class MergeVariantSitesIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        String fn = "merge1.vcf";

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
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.add(normalizePath(input));

        args.add("-V");
        File input2 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input2);
        args.add(normalizePath(input2));

        return args;
    }
}
