package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class Save10xBarcodesIntegrationTest extends BaseIntegrationTest {
    @Test
    public void basicTest() throws Exception{
        IntegrationTestSpec spec = new IntegrationTestSpec(
               " --input " + normalizePath(getInput()) +
                    " --cbOutput %s" +
                    " --umiOutput %s" +
                    " --tmp-dir " + getTmpDir(),
                Arrays.asList(
                    getTestFile("cbOutput.txt").getPath(),
                    getTestFile("umiOutput.txt").getPath())
                );

        spec.executeTest("basicTest", this);
    }

    private File getInput(){
        return new File(testBaseDir, "Population_iv.overlapping.bam");
    }
}
