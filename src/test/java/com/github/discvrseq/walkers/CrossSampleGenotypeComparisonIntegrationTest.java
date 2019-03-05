package com.github.discvrseq.walkers;

import htsjdk.tribble.bed.BEDCodec;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class CrossSampleGenotypeComparisonIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(
                    getTestFile("crossSampleGenotypeComparison.txt").getPath(),
                    getTestFile("crossSampleGenotypeComparisonSummary.txt").getPath()
                ));

        spec.executeTest("testBasicOperation", this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--variant");
        File input = new File(testBaseDir, "MendelianViolationEval.vcf");
        ensureVcfIndex(input);

        args.add(normalizePath(input));

        args.add("-referenceVCF");
        args.add(normalizePath(input));

        args.add("-O");
        args.add("%s");
        args.add("-summaryTable");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        return args;
    }
}
