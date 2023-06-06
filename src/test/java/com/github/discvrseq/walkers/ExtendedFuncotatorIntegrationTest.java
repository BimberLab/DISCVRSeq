package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.codecs.table.TableCodec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class ExtendedFuncotatorIntegrationTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "simpleTest.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(
                    getTestFile("funcotator1.vcf").getPath()
                ));

        spec.executeTest("testBasicOperation", this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("R", normalizePath(getHg19Micro()));
        args.add("ref-version", "hg19");

        args.add("data-sources-path", normalizePath(new File(testBaseDir, "funcotator")));
        args.add("output-file-format", "VCF");
        args.add("cf", normalizePath(new File(testBaseDir, "funcotator/funcotatorFields.txt")));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        ensureIndex(new File(testBaseDir, "funcotator/testSource/hg19/testSource.table"), new TableCodec());

        return args;
    }
}
