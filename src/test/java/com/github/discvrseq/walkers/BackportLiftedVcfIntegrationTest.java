package com.github.discvrseq.walkers;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Arrays;

public class BackportLiftedVcfIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/walkers/BackportLiftedVcf");

    @Test
    public void doBasicTest() throws Exception {
        File fasta = getHg19Micro();

        IntegrationTestSpec spec = new IntegrationTestSpec(
            " -R " + normalizePath(fasta) +
            " -V " + normalizePath(getInputVcf()) +
            " --targetFasta " + normalizePath(fasta) +
            " -O " + "%s" +
            " --tmp-dir " + getTmpDir(),
            Arrays.asList(testBaseDir + "/backportLiftedOutput.vcf"));

        spec.executeTest("doBasicTest", this);
    }

    @Test
    public void doBasicTestBcftools() throws Exception {
        File fasta = getHg19Micro();

        File tempVcf = IOUtils.createTempFile("backportLiftedBcftools", ".vcf");
        try (BufferedReader reader = IOUtil.openFileForBufferedUtf8Reading(getInputVcf()); PrintWriter writer = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(tempVcf))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.replaceAll(BackportLiftedVcf.ORIGINAL_CONTIG, BackportLiftedVcf.ORIGINAL_CONTIG_BCF);
                line = line.replaceAll(BackportLiftedVcf.ORIGINAL_START, BackportLiftedVcf.ORIGINAL_START_BCF);
                line = line.replaceAll(BackportLiftedVcf.ORIGINAL_ALLELES, BackportLiftedVcf.ORIGINAL_ALLELES_BCF);

                writer.println(line);
            }
        }

        try {
            IntegrationTestSpec spec = new IntegrationTestSpec(
                    " -R " + normalizePath(fasta) +
                            " -V " + normalizePath(tempVcf) +
                            " --targetFasta " + normalizePath(fasta) +
                            " -O " + "%s" +
                            " --tmp-dir " + getTmpDir(),
                    Arrays.asList(testBaseDir + "/backportLiftedOutput.vcf"));

            spec.executeTest("doBasicTest", this);
        }
        finally {
            tempVcf.delete();
        }
    }

    private File getInputVcf(){
        return new File(testBaseDir, "backportLifted.vcf");
    }
}
