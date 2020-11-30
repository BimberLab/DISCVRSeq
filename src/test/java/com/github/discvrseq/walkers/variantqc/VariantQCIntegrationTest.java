package com.github.discvrseq.walkers.variantqc;

import com.github.discvrseq.walkers.BaseIntegrationTest;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class VariantQCIntegrationTest extends BaseIntegrationTest {
    @Test
    public void testBasicOperation() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-L");
        args.addRaw("1");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperation", this);
        expected.delete();
    }

    @Test
    public void testMaxContigs() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testMaxContigs.html"));
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-maxContigs");
        args.addRaw("1");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testMaxContigs", this);
        expected.delete();
    }

    @Test
    public void testBasicOperationWithRawData() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-L");
        args.addRaw("1");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("-rd");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath(), getTestFile("testBasicOperation.json").getPath()));

        spec.executeTest("testBasicOperationWithRawData", this);
        expected.delete();
    }

    private ArgumentsBuilder getBasePedigreeArgs()
    {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "MendelianViolationEval.vcf");
        args.addRaw(normalizePath(input));
        ensureVcfIndex(input);

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-ped");
        args.addRaw(normalizePath(new File(testBaseDir, "MendelianViolationEval.ped")));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        return args;
    }

    @Test
    public void testPedigreeOutput() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testPedigreeOutput.html"));
        ArgumentsBuilder args = getBasePedigreeArgs();
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testPedigreeOutput", this);
        expected.delete();
    }

    @Test
    public void testPedigreeOutputWithValidation() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testPedigreeOutput.html"));
        ArgumentsBuilder args = getBasePedigreeArgs();
        args.addRaw("-pedValidationType");
        args.addRaw("SILENT");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testPedigreeOutputWithValidation", this);
        expected.delete();
    }

    @Test
    public void testExtendedReports() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testExtendedReports.html"));
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-L");
        args.addRaw("1");

        File extendedReports = getTestFile("extendedReports.txt");
        args.addRaw("-arf");
        args.addRaw(normalizePath(extendedReports));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testExtendedReports", this);
        expected.delete();
    }

    @Test
    public void testBasicOperationNoSamples() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperationNoSamples.html"));
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "clinvarV2.vcf");
        args.addRaw(normalizePath(input));

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-L");
        args.addRaw("1");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperationNoSamples", this);
        expected.delete();
    }

    private File generateCompleteOutput(File input) throws IOException {
        File out = new File(getTmpDir(), input.getName() + ".complete.html");
        try (BufferedReader reader = IOUtil.openFileForBufferedUtf8Reading(input); PrintWriter writer = new PrintWriter(IOUtil.openFileForBufferedUtf8Writing(out))) {
            HtmlGenerator.printStaticContent(writer);
            String line;
            while ((line = reader.readLine()) != null) {
                writer.println(line);
            }
        }

        return out;
    }
}
