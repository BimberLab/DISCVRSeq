package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class ClipOverlappingAlignmentsIntegrationTest extends BaseIntegrationTest{
    @Test
    public void doBasicTest() throws Exception{
        File fasta = getTestFile("SIVmac239.fasta");
        File bam = getTestFile("clipInput.sam");
        File bed = getTestFile("SIVmac239-overlap.bed");

        File outFile = IOUtils.createTempFile("clipOverlappingAlignmentsIntegrationTest", ".sam");
        File outReport = IOUtils.createTempFile("clipOverlappingAlignmentsIntegrationTest", ".txt");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -I " + normalizePath(bam) +
                        " --clipIntervals " + normalizePath(bed) +
                        " -rf " + normalizePath(outReport) +
                        " -O " + normalizePath(outFile) +
                        " --tmp-dir " + getTmpDir(),
                Collections.emptyList());

        spec.executeTest("doBasicTest", this);

        File expectedOutput = getTestFile("expectedOutput.sam");
        SamAssertionUtils.assertSamsEqual(outFile, expectedOutput);

        File expectedReport = getTestFile("expectedOutput.txt");
        IntegrationTestSpec.assertEqualTextFiles(outReport, expectedReport);

        outFile.delete();
        outReport.delete();
    }
}
