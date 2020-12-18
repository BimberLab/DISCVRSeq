package com.github.discvrseq.walkers.tagpcr;

import com.github.discvrseq.walkers.BaseIntegrationTest;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;

public class TagPcrSummaryIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testYamlParsing() throws Exception {
        Resource yml1 = new Resource("pENTR-PB511-Repro.yml", TagPcrSummary.class);
        File tmp1 = createTempFile("piggyBac", ".yml");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(tmp1)) {
            IOUtils.copy(yml1.getResourceContentsAsStream(), writer, StandardCharsets.UTF_8);
        }

        Resource yml2 = new Resource("Lentivirus.yml", TagPcrSummary.class);
        File tmp2 = createTempFile("lentivirus", ".yml");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(tmp2)) {
            IOUtils.copy(yml2.getResourceContentsAsStream(), writer, StandardCharsets.UTF_8);
        }

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(getHg19Micro()) +
                        " --insert-definition " + normalizePath(tmp1) +
                        " --insert-definition " + normalizePath(tmp2) +
                        " --validate-descriptor-only " +
                        " --tmp-dir " + getTmpDir(),
                Collections.emptyList());

        spec.executeTest("testYamlParsing", this);

        tmp1.delete();
        tmp2.delete();
    }

    @Test
    public void testYamlParsing2() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(getHg19Micro()) +
                        " --insert-name lentivirus" +
                        " --insert-name piggyBac" +
                        " --validate-descriptor-only " +
                        " --tmp-dir " + getTmpDir(),
                Collections.emptyList());

        spec.executeTest("testYamlParsing2", this);
    }

    @Test
    public void testShowDescriptors() throws Exception {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(getHg19Micro()) +
                        " --write-default-descriptors %s" +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(getTestFile("descriptors.txt").getPath()));

        spec.executeTest("testShowDescriptors", this);
    }

    @Test
    public void doBasicTest1() throws Exception {
        doTest(3);
    }

    @Test
    public void doBasicTest2() throws Exception {
        doTest(1);
    }

    private void doTest(int minAlign) throws Exception {
        String name = "BasicTest";
        File bam = new File(testBaseDir, "tagPcrTest.sam");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + getHg19Micro() +
                        " --bam " + normalizePath(bam) +
                        " --output-table %s " +
                        " --metrics-table %s " +
                        " -ma " + minAlign + " " +
                        " --reads-to-output 3 " +
                        " --insert-name piggybac " +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(getTestFile(name + "-" + minAlign + ".outputTable.txt").getPath(), getTestFile(name + "-" + minAlign + ".metrics.txt").getPath())
        );

        spec.executeTest(name, this);

    }

    //TODO: test primer design:
//        " --genbank-output sites.gb " +
//        " --primer-pair-table primerTable.txt " +
//        " --primer3-path primer3_core.exe " +
//        " --blastn-path blastn.exe " +
//        " --blast-db-path /Blast/FDC9E311-9D8F-1038-A30B-F8F3FC862593 " +
}
