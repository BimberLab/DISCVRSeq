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
        args.add("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.add(normalizePath(input));

        File fasta = getHg19Micro();
        args.add("-R");
        args.add(normalizePath(fasta));

        args.add("-L");
        args.add("1");

        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperation", this);
        expected.delete();
    }

    @Test
    public void testBasicOperationWithRawData() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--variant");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.add(normalizePath(input));

        File fasta = getHg19Micro();
        args.add("-R");
        args.add(normalizePath(fasta));

        args.add("-L");
        args.add("1");

        args.add("-O");
        args.add("%s");
        args.add("-rd");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath(), getTestFile("testBasicOperation.json").getPath()));

        spec.executeTest("testBasicOperationWithRawData", this);
        expected.delete();
    }

    private ArgumentsBuilder getBasePedigreeArgs()
    {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--variant");
        File input = new File(testBaseDir, "MendelianViolationEval.vcf");
        args.add(normalizePath(input));
        ensureVcfIndex(input);

        File fasta = getHg19Micro();
        args.add("-R");
        args.add(normalizePath(fasta));

        args.add("-ped");
        args.add(normalizePath(new File(testBaseDir, "MendelianViolationEval.ped")));

        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

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
        args.add("-pedValidationType");
        args.add("SILENT");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testPedigreeOutputWithValidation", this);
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
