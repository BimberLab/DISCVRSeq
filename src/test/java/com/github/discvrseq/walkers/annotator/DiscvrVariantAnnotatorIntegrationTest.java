package com.github.discvrseq.walkers.annotator;

import com.github.discvrseq.walkers.BaseIntegrationTest;
import org.broadinstitute.hellbender.exceptions.UserException;
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

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(normalizePath(getTestFile("/basicTestOutput.vcf"))));

        spec.executeTest("basicTest", this);
    }

    @Test
    public void basicTestWithCustomAnnotationsMissingArg() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File input = new File(testBaseDir, "MendelianViolationEval.vcf");
        args.add("V", normalizePath(input));

        File pedigree = new File(testBaseDir, "MendelianViolationEval.ped");
        args.add("ped", normalizePath(pedigree));

        args.add("A", "GenotypeConcordanceBySite");
        args.add("A", "MendelianViolationCount");
        args.add("A", "MinorAlleleFrequency");
        args.add("A", "Impact");

        args.add("A", "GenotypeConcordance");
        args.add("A", "MendelianViolationBySample");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), 1, UserException.BadInput.class);

        spec.executeTest("basicTestWithCustomAnnotationsMissingArg", this);
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

        args.add("A", "ChromosomeCounts");

        args.add("A", "GenotypeConcordance");
        args.add("A", "MendelianViolationBySample");

        args.addRaw("-rg");
        args.addRaw(normalizePath(input));

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

    @Test
    public void basicTestWithRefAlleleFrequency() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File input = new File(testBaseDir, "variantConcordanceRef1.vcf");
        ensureVcfIndex(input);
        args.add("V", normalizePath(input));

        args.add("A", "RefAlleleFrequency");
        args.add("target-info-field-key", "AF.2");

        args.addRaw("--af-source-vcf");
        File input2 = new File(testBaseDir, "variantConcordanceRef2.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(normalizePath(getTestFile("/basicTestWithRefAlleleFrequency.vcf"))));

        spec.executeTest("basicTestWithRefAlleleFrequency", this);
    }
}
