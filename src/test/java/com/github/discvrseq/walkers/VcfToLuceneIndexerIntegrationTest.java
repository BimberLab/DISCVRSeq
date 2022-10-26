package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class VcfToLuceneIndexerIntegrationTest extends BaseIntegrationTest {

    // Case where all IF are valid
    @Test
    public void doBasicTest() throws Exception {
        ensureVcfIndex(getInputVcf());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                  " -V " + normalizePath(getInputVcf()) +
                       " -O " + "%s" +
                       " -IF AD" +
                       " -IF PURPOSE" +
                       " -AN SampleList " +
                       " --tmp-dir " + getTmpDir(),
                Arrays.asList(normalizePath(getTestFile( "vcfFilterComparison.txt"))));

        spec.executeTest("doBasicTest", this);
    }

    // Case where there is an invalid IF
    @Test
    public void doInvalidIFTest() throws Exception {
        ensureVcfIndex(getInputVcf());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                  " -V " + normalizePath(getInputVcf()) +
                       " -O " + "%s" +
                       " -IF AD" +
                       " -IF PURPOSE" +
                       " -IF Foo" +  // does not exist in VCF
                       " -AN SampleList " +
                       " --tmp-dir " + getTmpDir(),
                       1,
                       GATKException.class);

        spec.executeTest("doInvalidIFTest", this);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "ClinvarAnnotator.vcf");
    }
}
