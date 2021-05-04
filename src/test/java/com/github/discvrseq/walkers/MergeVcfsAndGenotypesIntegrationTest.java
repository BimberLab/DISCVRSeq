package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class MergeVcfsAndGenotypesIntegrationTest extends BaseIntegrationTest {

    @Test
    public void basicTest() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("--variant:a");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-V:b");
        File input2 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-priority");
        args.addRaw("a,b");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        String fn = "basicTestOutput.vcf";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest("basicTest", this);
    }

    @Test
    public void basicTest2() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("--variant:a");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-V:b");
        File input2 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-V:c");
        File input3 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input3);
        args.addRaw(normalizePath(input3));

        args.add("genotypeMergeOption", "PRIORITIZE");

        args.addRaw("-priority");
        args.addRaw("a,b,c");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        String fn = "basicTest2Output.vcf";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest("basicTest", this);
    }

    @Test
    public void basicTest3() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("--variant:a");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-V:b");
        File input2 = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-priority");
        args.addRaw("a,b");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                1, UserException.BadInput.class);

        spec.executeTest("basicTest", this);


        args.add("genotypeMergeOption", "PRIORITIZE");

        String fn = "basicTest3Output.vcf";
        IntegrationTestSpec spec2 = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        executeTest(spec, "basicTest");
    }

    private void executeTest(IntegrationTestSpec spec, String name) throws IOException {
        try {
            spec.executeTest(name, this);
        }
        catch (AssertionError e) {
            throw e;
        }
    }
}
