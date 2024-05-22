package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class IdentifySoftClippedLociIntegrationTest extends BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        doTest("testBasicOperation", args, "IdentifySoftClippedLoci.bed");
    }

    @Test
    public void testWithFilter() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.addRaw("--minFraction");
        args.addRaw("0.2");

        doTest("testWithFilter", args, "IdentifySoftClippedLociFiltered.bed");
    }

    private void doTest(String name, ArgumentsBuilder args, String expectedFile) throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(expectedFile).getPath()));

        spec.executeTest(name, this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File fasta  = new File(testBaseDir, "Rhesus_KIR_and_MHC_1.0.fasta");
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-I");
        args.addRaw(normalizePath(new File(testBaseDir, "ImmunoGenotyper.qsort.bam")));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }

    @Test
    public void doTestWithVcf() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("-I");
        args.addRaw(normalizePath(new File(testBaseDir, "SoftClipTestData.sam")));

        File vcf = new File(testBaseDir, "softClipVariants.vcf");
        ensureVcfIndex(vcf);
        args.addRaw("-rv");
        args.addRaw(normalizePath(vcf));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile("IdentifySoftClippedLociWithVcf.bed").getPath()));

        spec.executeTest("doTestWithVcf", this);
    }
}
