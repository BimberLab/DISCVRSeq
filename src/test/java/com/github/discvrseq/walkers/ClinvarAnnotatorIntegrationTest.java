package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class ClinvarAnnotatorIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        doTest("testBasicOperation", args, "ClinvarAnnotatorOutput.vcf");
    }

    private void doTest(String name, ArgumentsBuilder args, String fn) throws Exception {
        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        File fasta = downloadHg19Micro();
        args.add("-R");
        args.add(normalizePath(fasta));

        IntegrationTestSpec spec = new IntegrationTestSpec(
            args.getString(),
            Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest(name, this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File testBaseDir = new File(publicTestDir + "com/github/discvrseq/TestData");

        //Must also create index once for clinvar VCF
        args.add("--clinvar");
        File clinvar = new File(testBaseDir, "clinvarV2.vcf");
        ensureVcfIndex(clinvar);
        args.add(normalizePath(clinvar));

        args.add("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.add(normalizePath(input));

        return args;
    }
}
