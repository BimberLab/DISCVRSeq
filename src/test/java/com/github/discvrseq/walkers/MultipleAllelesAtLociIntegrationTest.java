package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class MultipleAllelesAtLociIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        doTest("testBasicOperation", args, "MultipleAllelesAtLoci.bed");
    }

    @Test
    public void testWithFilter() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("--minBasePct");
        args.add("0.2");

        doTest("testWithFilter", args, "MultipleAllelesAtLociFiltered.bed");
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
        args.add("-R");
        args.add(normalizePath(fasta));

        args.add("-I");
        args.add(normalizePath(new File(testBaseDir, "ImmunoGenotyper.qsort.bam")));

        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        return args;
    }
}
