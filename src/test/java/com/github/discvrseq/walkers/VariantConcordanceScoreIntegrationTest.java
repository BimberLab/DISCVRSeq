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
        args.add("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.add(normalizePath(input));

        args.add("--ref-sites:REF1");
        File rs1 = new File(testBaseDir, "basicVcf.vcf");
        ensureVcfIndex(rs1);
        args.add(normalizePath(rs1));

        args.add("-rs:REF2");
        File rs2 = new File(testBaseDir, "clinvarV2.vcf");
        args.add(normalizePath(rs2));

        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        File expected = getTestFile("variantConcordanceOutput.txt");
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperation", this);
    }
}
