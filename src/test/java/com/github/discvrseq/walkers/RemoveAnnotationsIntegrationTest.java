package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class RemoveAnnotationsIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("-A");
        args.add("MAC");
        args.add("-XA");  //-A will trump this
        args.add("PURPOSE");
        args.add("-GA");
        args.add("DP");
        args.add("-XGA");
        args.add("AD");
        args.add("-ef");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile("testBasicOperation.vcf").getPath()));

        spec.executeTest("testBasicOperation", this);
    }

    @Test
    public void testBasicOperationSitesOnly() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("-A");
        args.add("MAC");
        args.add("-XA");
        args.add("PURPOSE");
        args.add("-GA");
        args.add("DP");
        args.add("-XGA");
        args.add("AD");
        args.add("-ef");
        args.add("-cgf");
        args.add("--sitesOnly");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile("testBasicOperationSitesOnly.vcf").getPath()));

        spec.executeTest("testBasicOperationSitesOnly", this);
    }

    @Test
    public void testExcludeSiteAnnotation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("-XA");
        args.add("PURPOSE");
        args.add("-GA");
        args.add("DP");
        args.add("-XGA");
        args.add("AD");
        args.add("-ef");
        args.add("-cgf");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile("testExcludeSite.vcf").getPath()));

        spec.executeTest("testExcludeSite", this);
    }


    private ArgumentsBuilder getBaseArgs()
    {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-R");
        args.add(normalizePath(getHg19Micro()));
        args.add("-V");
        args.add(normalizePath(new File(testBaseDir, "ClinvarAnnotator.vcf")));
        args.add("-O %s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        return args;
    }
}
