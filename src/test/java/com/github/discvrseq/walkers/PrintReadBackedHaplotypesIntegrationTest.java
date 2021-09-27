package com.github.discvrseq.walkers;

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class PrintReadBackedHaplotypesIntegrationTest extends BaseIntegrationTest {
    private void doExecute(IntegrationTestSpec spec, String name, CommandLineProgramTester test) throws IOException {
        try {
            spec.executeTest(name, test);
        }
        catch (Exception e) {
            spec.expectedFileNames();

            throw e;
        }
    }

    @Test
    public void filteredTest() throws Exception{
        File fasta = getGenome();

        File intervals = new File(getTmpDir(), "test.intervals");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(intervals))
        {
            writer.write("080_pRR_Reporter_U24:3566-3685\n");
            writer.write("080_pRR_Reporter_U24:3864-3983\n");
        }

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -I " + normalizePath(getInput()) +
                        " -O " + "%s" +
                        " -rc " + "0.75" +
                        " -mq " + "20" +
                        " -L " + normalizePath(intervals) +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(getTestFile("filteredTest.txt").getPath()));

        doExecute(spec, "filteredTest", this);

        intervals.delete();
    }

    @Test
    public void filteredTest2() throws Exception{
        File fasta = getGenome();

        File intervals = new File(getTmpDir(), "test.intervals");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(intervals))
        {
            writer.write("080_pRR_Reporter_U24:3566-3685\n");
            writer.write("080_pRR_Reporter_U24:3864-3983\n");
        }

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -I " + normalizePath(getInput()) +
                        " -O " + "%s" +
                        " -rc " + "0.95" +
                        " -mq " + "20" +
                        " -mmq " + "20" +
                        " -L " + normalizePath(intervals) +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(getTestFile("filteredTest2.txt").getPath()));

        doExecute(spec, "filteredTest2", this);

        intervals.delete();
    }

    @Test
    public void basicTest() throws Exception{
        File fasta = getGenome();

        File intervals = new File(getTmpDir(), "test.intervals");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(intervals))
        {
            writer.write("080_pRR_Reporter_U24:3566-3685\n");
            writer.write("080_pRR_Reporter_U24:3864-3983\n");
        }

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -I " + normalizePath(getInput()) +
                        " -O " + "%s" +
                        " -L " + normalizePath(intervals) +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(getTestFile("basicTest.txt").getPath()));

        doExecute(spec, "basicTest", this);

        intervals.delete();
    }

    @Test
    public void doTestWithoutIntervals() throws Exception{
        File fasta = getHg19Micro();

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -I " + normalizePath(getInput()) +
                        " -O " + "%s" +
                        " --tmp-dir " + getTmpDir(),
                1,
                CommandLineException.MissingArgument.class);

        spec.executeTest("doTestWithoutIntervals", this);
    }

    private File getGenome(){
        return new File(testBaseDir, "155_pRR_Reporter_U24.fasta");
    }

    private File getInput(){
        return new File(testBaseDir, "Population_vii.overlapping.bam");
    }
}
