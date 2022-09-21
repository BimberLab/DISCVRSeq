package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class VcfToLuceneIndexerIntegrationTest extends BaseIntegrationTest {
    @Test
    public void doBasicTest() throws Exception {
        ensureVcfIndex(getInputVcf());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                  " -V " + normalizePath(getInputVcf()) +
                       " -O " + "%s" +
                       " -IF AD" +
                       " -IF PURPOSE" +
                       " -IF Foo" +  // does not exist in VCF
                       " -AN SampleList " +
                       " --tmp-dir " + getTmpDir(),

                Arrays.asList(normalizePath(getTestFile( "vcfFilterComparison.txt"))));

        spec.executeTest("doBasicTest", this);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "ClinvarAnnotator.vcf");
    }
}
