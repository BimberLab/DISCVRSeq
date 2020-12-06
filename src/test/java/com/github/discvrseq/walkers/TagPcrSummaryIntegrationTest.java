package com.github.discvrseq.walkers;

import com.github.discvrseq.walkers.tagpcr.TagPcrSummary;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.charset.StandardCharsets;
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
        File tmp2 = createTempFile("defaults", ".txt");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(getHg19Micro()) +
                        " --write-default-descriptors " + normalizePath(tmp2) +
                        " --tmp-dir " + getTmpDir(),
                Collections.emptyList());

        spec.executeTest("testShowDescriptors", this);
    }
}
