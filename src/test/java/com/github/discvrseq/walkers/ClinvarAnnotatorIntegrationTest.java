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
        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        IntegrationTestSpec spec = new IntegrationTestSpec(
            args.getString(),
            Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest(name, this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        //Must also create index once for clinvar VCF
        args.addRaw("--clinvar");
        File clinvar = new File(testBaseDir, "clinvarV2.vcf");
        ensureVcfIndex(clinvar);
        args.addRaw(normalizePath(clinvar));

        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        return args;
    }
}
