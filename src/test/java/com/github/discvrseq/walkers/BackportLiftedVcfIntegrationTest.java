package com.github.discvrseq.walkers;

import com.github.discvrseq.TestUtils;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

public class BackportLiftedVcfIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/walkers/BackportLiftedVcf");

//    final ReferenceSequence reference = new ReferenceSequence("chr1", 0,
//            "CAAAAAAAAAACGTACGTACTCTCTCTCTACGT".getBytes());
//    //       123456789 123456789 123456789 123

    @Test
    public void doBasicTest() throws Exception{
        ArgumentsBuilder args = getBaseArgs();

        args.add("-O");
        File outFile = getSafeNonExistentFile("backportLiftedOutput.vcf");
        args.add(TestUtils.fixFilePath(outFile));

        runCommandLine(args.getArgsArray());

        File expected = new File(testBaseDir, "backportLiftedOutput.vcf");
        IntegrationTestSpec.assertEqualTextFiles(new File(TestUtils.fixFilePath(outFile)), expected);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "backportLifted.vcf");
    }

    private ArgumentsBuilder getBaseArgs() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File fasta = downloadHg19Micro();

        args.add("-R");
        args.add(TestUtils.fixFilePath(fasta));

        args.add("-V");
        args.add(TestUtils.fixFilePath(getInputVcf()));

        args.add("--targetFasta");
        args.add(TestUtils.fixFilePath(fasta));

        return args;
    }
}
