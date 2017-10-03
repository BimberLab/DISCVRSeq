package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

import static com.github.discvrseq.walkers.ImmunoGenotyper.*;

public class ImmunoGenotyperIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();
        args.add("--requireValidPair");

        doTest(args, "ImmunoGenotyperOutput");
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

        doTest(args, "ImmunoGenotyperOutputMM");
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

        doTest(args, "ImmunoGenotyperOutputNVP");
    }

    private void doTest(ArgumentsBuilder args, String fn) throws Exception{
        args.add("-O");

        File outFile = getSafeNonExistentFile(fn);
        String outFilePrefix = fixFilePath(outFile);
        args.add(outFilePrefix);

        runCommandLine(args.getArgsArray());

        for (String extention : Arrays.asList(GENOTYPE_EXTENSION, SUMMARY_EXTENSION, MISMATCH_EXTENSION)){
            File expected = getTestFile(fn + extention);
            File actual = IOUtils.getPath(outFilePrefix + extention).toFile();
            IntegrationTestSpec.assertEqualTextFiles(actual, expected);
        }
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File testBaseDir = new File(publicTestDir + "com/github/discvrseq/TestData");
        args.add("-R");
        args.add(fixFilePath(new File(testBaseDir, "Rhesus_KIR_and_MHC_1.0.fasta")));

        args.add("-I");
        args.add(fixFilePath(new File(testBaseDir, "ImmunoGenotyper.qsort.bam")));

        args.add("--referenceToLineageFile");
        args.add(fixFilePath(new File(testBaseDir, "lineageMap.txt")));

        return args;
    }
}
