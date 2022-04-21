package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class AnnotateNovelSitesTest extends  BaseIntegrationTest {

    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = getBaseArgs();

        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-rv");
        File ref = new File(testBaseDir, "AnnotateNovelSitesRef.vcf");
        ensureVcfIndex(ref);
        args.addRaw(normalizePath(ref));

        args.add("an", "mGAPV");
        args.add("ad", "'The first mGAP version where variants at this site appeared'");
        args.add("av", "2.0");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(
                    getTestFile("output.vcf").getPath(),
                    getTestFile("missing.vcf").getPath(),
                        getTestFile("novel.vcf").getPath()
                ));

        spec.executeTest("testBasicOperation", this);
    }

    private ArgumentsBuilder getBaseArgs() {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("-ns");
        args.addRaw("%s");
        args.addRaw("-ms");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }
}
