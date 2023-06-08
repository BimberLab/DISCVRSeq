package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class SampleSpecificGenotypeFiltrationIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs("ClinvarAnnotator.vcf");

        String fn = "basicTestOut.vcf";

        args.addFlag(VariantFiltration.NO_CALL_GTS_LONG_NAME);
        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        args.add(VariantFiltration.GENOTYPE_FILTER_EXPRESSION_LONG_NAME, "Set1:DP<15");
        args.add(VariantFiltration.GENOTYPE_FILTER_NAME_LONG_NAME, "DP-LT15");

        args.add(VariantFiltration.GENOTYPE_FILTER_EXPRESSION_LONG_NAME, "Set2:DP<10");
        args.add(VariantFiltration.GENOTYPE_FILTER_NAME_LONG_NAME, "DP-LT10");

        args.add(VariantFiltration.GENOTYPE_FILTER_EXPRESSION_LONG_NAME, "Set2:DP>30");
        args.add(VariantFiltration.GENOTYPE_FILTER_NAME_LONG_NAME, "DP-GT30");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Collections.singletonList(getTestFile(fn).getPath()));

        spec.executeTest("testBasicOperation", this);
    }

    @Test
    public void testRGQFilter() throws Exception {
        ArgumentsBuilder args = getBaseArgs("genotypeFilterTest.vcf");

        String fn = "rgqTestOut.vcf";

        args.addFlag(VariantFiltration.NO_CALL_GTS_LONG_NAME);
        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        args.add(VariantFiltration.GENOTYPE_FILTER_EXPRESSION_LONG_NAME, "'Set2:GQ<20 && RGQ<20'");
        args.add(VariantFiltration.GENOTYPE_FILTER_NAME_LONG_NAME, "GQ-LT20");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Collections.singletonList(getTestFile(fn).getPath()));

        spec.executeTest("testBasicOperation", this);
    }

    private ArgumentsBuilder getBaseArgs(String inputVcf) {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));
        args.addRaw("--variant");
        File input = new File(testBaseDir, inputVcf);
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.add("sample-map", normalizePath(new File(testBaseDir, "sample-map.txt")));
        return args;
    }
}
