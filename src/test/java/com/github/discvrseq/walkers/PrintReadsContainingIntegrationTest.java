package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PrintReadsContainingIntegrationTest extends BaseIntegrationTest {
    @DataProvider(name = "testExpressionsData")
    public Object[][] testExpressionsData() {
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{"test1", new String[]{"TGGTGAAACCCTGTCTCT"}, null, null, false, true});
        tests.add(new Object[]{"test2", new String[]{"^CTTATCCTGTGGCTGCTTGA$"}, null, null, false, true});
        tests.add(new Object[]{"test2", null, new String[]{"^CTTATCCTGTGGCTGCTTGA$"}, null, false, true});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "testExpressionsData")
    public void testExpressionsPaired(String testName, String[] exprs, String[] r1Exprs, String[] r2Exprs, boolean matchAllExpressions, boolean paired) throws IOException {

        ArgumentsBuilder args = getBaseArgs(paired);

        if (exprs != null) {
            Arrays.stream(exprs).forEach(x -> {
                args.add("-e");
                args.add(x);
            });
        }

        if (r1Exprs != null) {
            Arrays.stream(r1Exprs).forEach(x -> {
                args.add("-e1");
                args.add(x);
            });
        }

        if (r2Exprs != null) {
            Arrays.stream(r2Exprs).forEach(x -> {
                args.add("-e2");
                args.add(x);
            });
        }

        if (matchAllExpressions) {
            args.add("-ma");
        }

        List<String> expected = new ArrayList<>();
        expected.add(getTestFile(testName + "_R1.fastq").getPath());
        if (paired) {
            expected.add(getTestFile(testName + "_R2.fastq").getPath());
        }

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                expected);

        try {
            spec.executeTest(testName, this);
        }
        catch (AssertionError e) {
            throw e;
        }
    }

    private ArgumentsBuilder getBaseArgs(boolean paired) {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--fastq");
        args.add(normalizePath(getTestFile("fq1.fastq")));

        args.add("--output");
        args.add("%s");

        if (paired) {
            args.add("--fastq2");
            args.add(normalizePath(getTestFile("fq2.fastq")));

            args.add("--output2");
            args.add("%s");
        }

        args.add("--tmp-dir");
        args.add(getTmpDir());




        return args;
    }
}
