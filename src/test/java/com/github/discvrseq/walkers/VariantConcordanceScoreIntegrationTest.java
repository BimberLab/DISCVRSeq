package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class VariantConcordanceScoreIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        args.addRaw("--ref-sites:REF1");
        File rs1 = new File(testBaseDir, "basicVcf.vcf");
        ensureVcfIndex(rs1);
        args.addRaw(normalizePath(rs1));

        args.addRaw("-rs:REF2");
        File rs2 = new File(testBaseDir, "clinvarV2.vcf");
        args.addRaw(normalizePath(rs2));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        File expected = getTestFile("variantConcordanceOutput.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperation", this);
    }
}
