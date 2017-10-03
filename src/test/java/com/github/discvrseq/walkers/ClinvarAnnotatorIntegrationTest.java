package com.github.discvrseq.walkers;

import com.github.discvrseq.Main;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class ClinvarAnnotatorIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        doTest(args, "ClinvarAnnotatorOutput.vcf");
    }

    private void doTest(ArgumentsBuilder args, String fn) throws Exception{
        args.add("-O");
        File outFile = getSafeNonExistentFile(fn);
        args.add(outFile);

        runCommandLine(args.getArgsArray());

        File expected = getTestFile(fn);
        IntegrationTestSpec.assertEqualTextFiles(outFile, expected);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File testBaseDir = new File(publicTestDir + "com/github/discvrseq/TestData");

        //Must also create index once for clinvar VCF
        args.add("--clinvar");
        File clinvar = new File(testBaseDir, "clinvarV2.vcf");
        args.add(clinvar);

        args.add("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.add(input);

        return args;
    }

    @Override
    public Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }
}
