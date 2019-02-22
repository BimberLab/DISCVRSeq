package com.github.discvrseq.walkers;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class PrintReadBackedHaplotypesIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/walkers/PrintReadBackedHaplotypes");

//    @Test
//    public void doBasicTest() throws Exception{
//        //TODO: rework this test
//        File fasta = new File(testBaseDir, "80_Mmul_8.0.1.fasta");
//
//        //chr14:69328275-69328393
//        IntegrationTestSpec spec = new IntegrationTestSpec(
//                " -R " + normalizePath(fasta) +
//                        " -I " + normalizePath(getInput()) +
//                        " -O " + "%s" +
//                        " --tmp-dir " + getTmpDir() +
//                        " -L " + normalizePath(new File(testBaseDir, "HaploTest.intervals")),
//                Arrays.asList(testBaseDir + "/haplotypeOutput.txt"));
//
//        spec.executeTest("doBasicTest", this);
//    }
//
//    @Test
//    public void doBasicTestFiltered() throws Exception {
//        File fasta = new File(testBaseDir, "80_Mmul_8.0.1.fasta");
//
//        //chr14:69328275-69328393
//        IntegrationTestSpec spec = new IntegrationTestSpec(
//                " -R " + normalizePath(fasta) +
//                        " -I " + normalizePath(getInput()) +
//                        " -O " + "%s" +
//                        " --tmp-dir " + getTmpDir() +
//                        " -rc 0.8" +
//                        " -L " + normalizePath(new File(testBaseDir, "HaploTest.intervals")),
//                Arrays.asList(testBaseDir + "/haplotypeOutputFiltered.txt"));
//
//        doExecute(spec, "doBasicTestFiltered", this);
//    }
//
//    @Test
//    public void doBasicTestFilteredBWA() throws Exception {
//        File fasta = new File(testBaseDir, "80_Mmul_8.0.1.fasta");
//
//        //chr14:69328275-69328393
//        File input = new File(testBaseDir, "WGA436-BWAmem.bam");
//        IntegrationTestSpec spec = new IntegrationTestSpec(
//                " -R " + normalizePath(fasta) +
//                        " -I " + normalizePath(input) +
//                        " -O " + "%s" +
//                        " --tmp-dir " + getTmpDir() +
//                        " -rc 0.8" +
//                        " -L " + normalizePath(new File(testBaseDir, "HaploTest.intervals")),
//                Arrays.asList(testBaseDir + "/haplotypeOutputFilteredBWA.txt"));
//
//        doExecute(spec, "doBasicTestFilteredBWA", this);
//    }

    private void doExecute(IntegrationTestSpec spec, String name, CommandLineProgramTester test) throws IOException {
        try {
            spec.executeTest(name, test);
        }
        catch (AssertionError e) {
            spec.expectedFileNames();

            throw e;
        }
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

    private File getInput(){
        return new File(testBaseDir, "WGA436.bam");
    }
}
