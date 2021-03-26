package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class Summarize10xAlignmentsIntegrationTest extends BaseIntegrationTest {

    @Test
    public void doBasicTest() throws Exception {
        File bam = new File(testBaseDir, "10xTestData.sam");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -I " + normalizePath(bam) +
                        " -O " + "%s" +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(normalizePath(getTestFile( "basicTest.txt"))));

        spec.executeTest("doBasicTest", this);
    }
}
