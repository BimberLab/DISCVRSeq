package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class BackportLiftedVcfIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/walkers/BackportLiftedVcf");

    @Test
    public void doBasicTest() throws Exception{
        File fasta = downloadHg19Micro();

        IntegrationTestSpec spec = new IntegrationTestSpec(
            " -R " + normalizePath(fasta) +
            " -V " + normalizePath(getInputVcf()) +
            " --targetFasta " + normalizePath(fasta) +
            " -O " + "%s" +
            " --tmp-dir " + getTmpDir(),
            Arrays.asList(testBaseDir + "/backportLiftedOutput.vcf"));

        spec.executeTest("doBasicTest", this);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "backportLifted.vcf");
    }
}
