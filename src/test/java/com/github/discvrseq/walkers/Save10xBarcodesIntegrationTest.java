package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class Save10xBarcodesIntegrationTest extends BaseIntegrationTest {
    @Test
    public void basicTest() throws Exception{
        IntegrationTestSpec spec = new IntegrationTestSpec(
               " --input " + normalizePath(getInput()) +
                    " --output %s" +
                    " --tmp-dir " + getTmpDir(),
                Collections.singletonList(
                        getTestFile("output.txt").getPath()
                ));

        spec.executeTest("basicTest", this);
    }

    private File getInput(){
        return new File(testBaseDir, "Population_iv.overlapping.bam");
    }
}
