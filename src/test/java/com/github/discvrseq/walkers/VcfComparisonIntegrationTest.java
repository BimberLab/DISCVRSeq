package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class VcfComparisonIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-referenceVCF");
        File ref = new File(testBaseDir, "clinvarV2.vcf");
        ensureVcfIndex(ref);
        args.addRaw(normalizePath(ref));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(
                    getTestFile("vcfComparison.txt").getPath(),
                    getTestFile("vcfComparisonSites.txt").getPath(),
                    getTestFile("novelSites.vcf").getPath(),
                    getTestFile("missingSites.vcf").getPath()
                ));

        spec.executeTest("testBasicOperation", this);
    }

    @Test
    public void testBasicOperationWithIncludeNoCall() throws Exception {
        ArgumentsBuilder args = getBaseArgs(false);

        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("--include-no-call-sites");

        args.addRaw("-referenceVCF");
        File ref = new File(testBaseDir, "clinvarV2.vcf");
        ensureVcfIndex(ref);
        args.addRaw(normalizePath(ref));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(
                        getTestFile("vcfComparisonWithNoCall.txt").getPath()
                ));

        spec.executeTest("testBasicOperation", this);
    }

    @Test
    public void testBasicOperation2() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-referenceVCF");
        File ref = new File(testBaseDir, "mergeVcf3.vcf");
        ensureVcfIndex(ref);
        args.addRaw(normalizePath(ref));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(
                        getTestFile("vcfComparison2.txt").getPath(),
                        getTestFile("vcfComparisonSites2.txt").getPath(),
                        getTestFile("novelSites2.vcf").getPath(),
                        getTestFile("missingSites2.vcf").getPath()
                ));

        spec.executeTest("testBasicOperation2", this);
    }

    @Test
    public void testFailure() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-referenceVCF");
        File ref = new File(testBaseDir, "clinvarV2.vcf");
        ensureVcfIndex(ref);
        args.addRaw(normalizePath(ref));

        args.add("V", normalizePath(new File(testBaseDir, "clinvarV2.vcf")));

        IntegrationTestSpec spec = new IntegrationTestSpec(args.getString(), 4, UserException.BadInput.class);

        spec.executeTest("testFailure", this);
    }

    private ArgumentsBuilder getBaseArgs() {
        return getBaseArgs(true);
    }

    private ArgumentsBuilder getBaseArgs(boolean includeFileOutputs) {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("-O");
        args.addRaw("%s");

        if (includeFileOutputs) {
            args.addRaw("-sites-output");
            args.addRaw("%s");

            args.addRaw("--novel-or-altered-sites-vcf");
            args.addRaw("%s");

            args.addRaw("--missing-sites-vcf");
            args.addRaw("%s");
        }

        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }
}
