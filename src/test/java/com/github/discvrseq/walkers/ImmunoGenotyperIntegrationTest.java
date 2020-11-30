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
        args.addRaw("--requireValidPair");

        doTest("testBasicOperation", args, "ImmunoGenotyperOutput");
    }

    @Test
    public void testWithMismatches() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.addRaw("--mismatchesTolerated");
        args.addRaw(10);
        args.addRaw("-mmq");
        args.addRaw(0);
        args.addRaw("--minReadCountForRef");
        args.addRaw(1);
        args.addRaw("--minPctForExport");
        args.addRaw(0.005);
        args.addRaw("--minReadCountForExport");
        args.addRaw("1");

        doTest("testWithMismatches", args, "ImmunoGenotyperOutputMM");
    }

    @Test
    public void testWithoutRequireValidPair() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.addRaw("-mmq");
        args.addRaw(0);
        args.addRaw("--minPctForRef");
        args.addRaw(0.001);

        args.addRaw("--minPctForExport");
        args.addRaw(0.001);

        doTest("testWithoutRequireValidPair", args, "ImmunoGenotyperOutputNVP");
    }

    private void doTest(String name, ArgumentsBuilder args, String fn) throws Exception{
        System.setProperty("java.io.tmpdir", getTmpDir());  //windows hack
        File outFile = new File(normalizePath(getSafeNonExistentFile(fn)));
        String outFilePrefix = normalizePath(outFile);
        args.addRaw("-O");
        args.addRaw(outFilePrefix);

        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

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
        args.addRaw("-R");

        File fasta  = new File(testBaseDir, "Rhesus_KIR_and_MHC_1.0.fasta");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-I");
        args.addRaw(normalizePath(new File(testBaseDir, "ImmunoGenotyper.qsort.bam")));

        args.addRaw("--referenceToLineageFile");
        args.addRaw(normalizePath(new File(testBaseDir, "lineageMap.txt")));

        return args;
    }
}
