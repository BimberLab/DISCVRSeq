package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class DiscvrVariantAnnotatorIntegrationTest extends BaseIntegrationTest {
    @Test
    public void basicTest() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File input = new File(testBaseDir, "MendelianViolationEval.vcf");
        args.add("V", normalizePath(input));

        File pedigree = new File(testBaseDir, "MendelianViolationEval.ped");
        args.add("ped", normalizePath(pedigree));
        
        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(normalizePath(getTestFile("/basicTestOutput.vcf"))));

        try {
            spec.executeTest("basicTest", this);
        }
        catch (AssertionError e) {
            throw e;
        }
    }

    @Test
    public void basicTestWithCustomAnnotations() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File input = new File(testBaseDir, "MendelianViolationEval.vcf");
        args.add("V", normalizePath(input));

        File pedigree = new File(testBaseDir, "MendelianViolationEval.ped");
        args.add("ped", normalizePath(pedigree));

        args.add("A", "GenotypeConcordanceBySite");
        args.add("A", "MendelianViolationCount");
        args.add("A", "MinorAlleleFrequency");

        args.add("A", "GenotypeConcordance");
        args.add("A", "MendelianViolationBySample");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(normalizePath(getTestFile("/basicTestWithCustomAnnotationsOutput.vcf"))));

        try {
            spec.executeTest("basicTestWithCustomAnnotations", this);
        }
        catch (AssertionError e) {
            throw e;
        }
    }
}
