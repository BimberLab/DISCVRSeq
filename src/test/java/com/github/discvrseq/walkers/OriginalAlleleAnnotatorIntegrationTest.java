package com.github.discvrseq.walkers;

import com.github.discvrseq.TestUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

public class OriginalAlleleAnnotatorIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/");

    @Test
    public void doBasicTest() throws Exception{
        ArgumentsBuilder args = getBaseArgs();

        args.add("-O");
        File outFile = getSafeNonExistentFile("originalAlleleAnnotatorOutput.vcf");
        args.add(TestUtils.fixFilePath(outFile));

        runCommandLine(args.getArgsArray());

        File expected = new File(testBaseDir, "walkers/OriginalAlleleAnnotator/originalAlleleAnnotatorOutput.vcf");
        IntegrationTestSpec.assertEqualTextFiles(new File(TestUtils.fixFilePath(outFile)), expected);
    }

    @Test
    public void doOriginalAlleleTest() throws Exception{
        ArgumentsBuilder args = getBaseArgs();

        args.add("-O");
        File outFile = getSafeNonExistentFile("originalAlleleAnnotatorOutput.vcf");
        args.add(TestUtils.fixFilePath(outFile));

        runCommandLine(args.getArgsArray());

        File expected = new File(testBaseDir, "walkers/OriginalAlleleAnnotator/originalAlleleAnnotatorOutput.vcf");
        IntegrationTestSpec.assertEqualTextFiles(new File(TestUtils.fixFilePath(outFile)), expected);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "TestData/basicVcf.vcf");
    }

    private ArgumentsBuilder getBaseArgs() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        File fasta = downloadHg19Micro();

        args.add("-R");
        args.add(TestUtils.fixFilePath(fasta));

        args.add("-V");
        args.add(TestUtils.fixFilePath(getInputVcf()));

        return args;
    }
}
