package com.github.discvrseq.walkers;

import htsjdk.tribble.bed.BEDCodec;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class AppendGenotypesIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        String fn = "append1.vcf";

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest("testBasicOperation", this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "basicVcf.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-g");
        File geno1 = new File(testBaseDir, "genotypeAppend1.vcf");
        ensureVcfIndex(geno1);
        args.addRaw(normalizePath(geno1));

        args.addRaw("-g");
        File geno2 = new File(testBaseDir, "genotypeAppend2.vcf");
        ensureVcfIndex(geno2);
        args.addRaw(normalizePath(geno2));

        return args;
    }
}
