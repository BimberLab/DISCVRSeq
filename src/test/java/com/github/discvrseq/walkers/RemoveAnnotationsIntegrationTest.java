package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class RemoveAnnotationsIntegrationTest extends BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        File cassandra = new File(testBaseDir, "ClinvarAnnotator.vcf");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                    " -V " + normalizePath(cassandra) +
                    " -XA PURPOSE" +
                    " -XA RN" +
                    " -GA DP" +
                    " -GA AD" +
                    " -GA GQ" +
                    " -GA GT" +
                    " -GA PL" +
                    " -ef" +
                    " -cgf" +
                    " -O " + "%s" +
                    " --tmp-dir " + getTmpDir(),
            Arrays.asList(normalizePath(getTestFile("/testBasicOperation.vcf"))));

        spec.executeTest("testBasicOperation", this);
    }

    @Test
    public void testBasicOperationSitesOnly() throws Exception {
        File cassandra = new File(testBaseDir, "ClinvarAnnotator.vcf");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -V " + normalizePath(cassandra) +
                        " -XA PURPOSE" +
                        " -XA RN" +
                        " -GA DP" +
                        " -GA AD" +
                        " -GA GQ" +
                        " -GA GT" +
                        " -GA PL" +
                        " -ef" +
                        " -cgf" +
                        " --sitesOnly" +
                        " -O " + "%s" +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(normalizePath(getTestFile("/testBasicOperationSitesOnly.vcf"))));

        spec.executeTest("testBasicOperationSitesOnly", this);
    }
}
