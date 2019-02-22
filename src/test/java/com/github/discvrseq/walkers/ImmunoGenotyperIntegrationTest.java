package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

import static com.github.discvrseq.walkers.ImmunoGenotyper.*;


public class ImmunoGenotyperIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("--requireValidPair");

        doTest("testBasicOperation", args, "ImmunoGenotyperOutput");
    }

    @Test
    public void testWithMismatches() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("--mismatchesTolerated");
        args.add(10);
        args.add("-mmq");
        args.add(0);
        args.add("--minReadCountForRef");
        args.add(1);
        args.add("--minPctForExport");
        args.add(0.005);
        args.add("--minReadCountForExport");
        args.add("1");

        doTest("testWithMismatches", args, "ImmunoGenotyperOutputMM");
    }

    @Test
    public void testWithoutRequireValidPair() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("-mmq");
        args.add(0);
        args.add("--minPctForRef");
        args.add(0.001);

        args.add("--minPctForExport");
        args.add(0.001);

        doTest("testWithoutRequireValidPair", args, "ImmunoGenotyperOutputNVP");
    }

    private void doTest(String name, ArgumentsBuilder args, String fn) throws Exception{
        System.setProperty("java.io.tmpdir", getTmpDir());  //windows hack
        File outFile = new File(normalizePath(getSafeNonExistentFile(fn)));
        String outFilePrefix = normalizePath(outFile);
        args.add("-O");
        args.add(outFilePrefix);

        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Collections.emptyList());

        spec.executeTest(name, this);

        for (String extention : Arrays.asList(GENOTYPE_EXTENSION, SUMMARY_EXTENSION, MISMATCH_EXTENSION)){
            File expected = getTestFile(fn + extention);
            File actual = IOUtils.getPath(outFilePrefix + extention).toFile();
            IntegrationTestSpec.assertEqualTextFiles(actual, expected);
        }
    }

    private ArgumentsBuilder getBaseArgs() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-R");

        File fasta  = new File(testBaseDir, "Rhesus_KIR_and_MHC_1.0.fasta");
        args.add(normalizePath(fasta));

        args.add("-I");
        args.add(normalizePath(new File(testBaseDir, "ImmunoGenotyper.qsort.bam")));

        args.add("--referenceToLineageFile");
        args.add(normalizePath(new File(testBaseDir, "lineageMap.txt")));

        return args;
    }
}
