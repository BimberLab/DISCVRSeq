package com.github.discvrseq.walkers.variantqc;

import com.github.discvrseq.walkers.BaseIntegrationTest;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;

public class VariantQCIntegrationTest extends BaseIntegrationTest {
    private ArgumentsBuilder getBaseArgs(boolean limitToChr1) {
        return(getBaseArgs(limitToChr1, null, 1));
    }

    private ArgumentsBuilder getBaseArgs(boolean limitToChr1, String vcfName, int threads) {
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant" + (vcfName == null ? "" : ":" + vcfName));
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        if (limitToChr1) {
            args.addRaw("-L");
            args.addRaw("1");
        }

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        if (threads > 1) {
            args.add("threads", threads);
        }

        return args;
    }

    @Test
    public void testBasicOperation() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
        ArgumentsBuilder args = getBaseArgs(true);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperation", this);
        expected.delete();
    }

    @Test
    public void testBasicOperationThreaded() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
        ArgumentsBuilder args = getBaseArgs(true, null, 2);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testBasicOperationThreaded", this);
        expected.delete();
    }

    @Test
    public void testMultiVcf() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testMultiVcf.html"));
        ArgumentsBuilder args = getBaseArgs(true, "vcf1", 1);
        args.addRaw("--variant:vcf2");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testMultiVcf", this);
        expected.delete();
    }

    @Test
    public void testMultiVcfThreaded() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testMultiVcf.html"));
        ArgumentsBuilder args = getBaseArgs(true, "vcf1", 2);
        args.addRaw("--variant:vcf2");
        File input = new File(testBaseDir, "ClinvarAnnotator.vcf");
        args.addRaw(normalizePath(input));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testMultiVcfThreaded", this);
        expected.delete();
    }

    @Test
    public void testMaxContigs() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testMaxContigs.html"));
        ArgumentsBuilder args = getBaseArgs(false);

        args.addRaw("-maxContigs");
        args.addRaw("1");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testMaxContigs", this);
        expected.delete();
    }

    @Test
    public void testBasicOperationWithRawData() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
        ArgumentsBuilder args = getBaseArgs(true);

        args.addRaw("-rd");
        args.addRaw("%s");

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
        ArgumentsBuilder args = getBaseArgs(true);

        File extendedReports = getTestFile("extendedReports.txt");
        args.addRaw("-arf");
        args.addRaw(normalizePath(extendedReports));

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

    @DataProvider(name = "testPedigreeValidationData")
    private Object[][] testPedigreeValidationData() {
        return new Object[][] {
                new Object[]{PedigreeValidationType.SILENT, false},
                new Object[]{PedigreeValidationType.STRICT, true},
                new Object[]{null, true}
        };
    }

    @Test(dataProvider = "testPedigreeValidationData")
    public void testPedigreeValidation(PedigreeValidationType pvt, boolean expectFail) throws IOException {
        String pedFile = "MendelianViolationEval.ped";
        String name = "testPedigreeValidation";

        ArgumentsBuilder args = getBaseArgs(true);
        args.add("ped", normalizePath(new File(testBaseDir, pedFile)));

        if (pvt != null) {
            args.add("pedValidationType", pvt.name());
        }

        IntegrationTestSpec spec;
        if (expectFail) {
            spec = new IntegrationTestSpec(args.getString(), 1, UserException.class);
        }
        else {
            File expected = generateCompleteOutput(getTestFile("testBasicOperation.html"));
            spec = new IntegrationTestSpec(args.getString(), Collections.singletonList(expected.getPath()));
        }

        spec.executeTest(name, this);
    }

    @Test
    public void testMaxContigs2() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testMaxContigs2.html"));

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "vcfTwoContigs.vcf");
        args.addRaw(normalizePath(input));
        ensureVcfIndex(input);

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        args.add("maxContigs", 1);

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testMaxContigs2", this);
        expected.delete();
    }

    @Test
    public void testMaxContigsWithRecovery() throws Exception {
        File expected = generateCompleteOutput(getTestFile("testMaxContigs3.html"));

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--variant");
        File input = new File(testBaseDir, "vcfTwoContigs.vcf");
        args.addRaw(normalizePath(input));
        ensureVcfIndex(input);

        File fasta = getHg19Micro();
        args.addRaw("-R");
        args.addRaw(normalizePath(fasta));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        args.add("maxContigs", 1);
        args.add("contigsToRetain", "2");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(), Arrays.asList(expected.getPath()));

        spec.executeTest("testMaxContigs3", this);
        expected.delete();
    }
}
