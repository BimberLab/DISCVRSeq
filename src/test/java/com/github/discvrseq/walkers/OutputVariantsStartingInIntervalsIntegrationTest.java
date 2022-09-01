package com.github.discvrseq.walkers;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;

public class OutputVariantsStartingInIntervalsIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));
        args.addRaw("--variant");
        File input = new File(testBaseDir, "basicVcf.vcf");
        ensureVcfIndex(input);
        args.addRaw(normalizePath(input));

        args.addRaw("-L");
        File intervalFile = new File(getTmpDir(), "intervals.list");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(intervalFile)) {
            writer.write("1:62-100");
        }
        args.addRaw(normalizePath(intervalFile));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Collections.singletonList(getTestFile("basicOutput.vcf").getPath()));

        spec.executeTest("testBasicOperation", this);

        intervalFile.delete();
    }
}
