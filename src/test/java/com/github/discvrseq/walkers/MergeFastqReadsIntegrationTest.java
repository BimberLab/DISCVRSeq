package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class MergeFastqReadsIntegrationTest extends BaseIntegrationTest {
    @Test
    public void basicMergeTest() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-fq1");
        File fq1 = new File(new File(getToolTestDataDir()).getParentFile(), "PrintReadsContaining/fq1.fastq");
        args.add(normalizePath(fq1));

        args.add("-fq2");
        File fq2 = new File(new File(getToolTestDataDir()).getParentFile(), "PrintReadsContaining/fq2.fastq");
        args.add(normalizePath(fq2));

        args.add("-O");
        args.add("%s");

        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile("basicMergeTest.fastq").getPath()));

        spec.executeTest("basicMergeTest", this);
    }

    @Test
    public void mergeTestWithMinLength() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-fq1");
        File fq1 = new File(new File(getToolTestDataDir()).getParentFile(), "PrintReadsContaining/fq1.fastq");
        args.add(normalizePath(fq1));

        args.add("-fq2");
        File fq2 = new File(new File(getToolTestDataDir()).getParentFile(), "PrintReadsContaining/fq2.fastq");
        args.add(normalizePath(fq2));

        args.add("-O");
        args.add("%s");

        args.add("--minLength");
        args.add("100");

        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile("mergeTestWithMinLength.fastq").getPath()));

        spec.executeTest("mergeTestWithMinLength", this);
    }
}
