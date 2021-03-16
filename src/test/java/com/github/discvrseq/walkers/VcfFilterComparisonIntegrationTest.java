package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class VcfFilterComparisonIntegrationTest extends BaseIntegrationTest {
    @Test
    public void doBasicTest() throws Exception {
        File fasta = getHg19Micro();
        ensureVcfIndex(getInputVcf());

        IntegrationTestSpec spec = new IntegrationTestSpec(
            " -R " + normalizePath(fasta) +
            " -V:vcf1 " + normalizePath(getInputVcf()) +
            " -V:vcf2 " + normalizePath(getInputVcf()) +
            " -O " + "%s" +
            " --tmp-dir " + getTmpDir(),
            Arrays.asList(normalizePath(getTestFile( "vcfFilterComparison.txt"))));

        spec.executeTest("doBasicTest", this);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "basicVcfFiltered.vcf");
    }
}
