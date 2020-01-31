package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PrintReadsContainingIntegrationTest extends BaseIntegrationTest {
    @DataProvider(name = "testExpressionsData")
    public Object[][] testExpressionsData() {
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{"test1", new String[]{"TGGTGAAACCCTGTCTCT"}, null, null, false, 4, 0, null});
        tests.add(new Object[]{"test1", null, new String[]{"TGGTGAAACCCTGTCTCT"}, null, false, 0, 0, null});
        tests.add(new Object[]{"test1", null, null, new String[]{"TGGTGAAACCCTGTCTCT"}, false, 4, 0, null});

        tests.add(new Object[]{"test2", new String[]{"^CTTATCCTGTGGCTGCTTGA$"}, null, null, false, 4, 4, null});
        tests.add(new Object[]{"test2", null, new String[]{"^CTTATCCTGTGGCTGCTTGA$"}, null, false, 4, 4, null});

        tests.add(new Object[]{"test3", new String[]{"TCATACTCGGAGGAGCTGG", "ACCTCCACCTCCCGG"}, null, null, false, 4, 4, 3});
        tests.add(new Object[]{"test3", null, new String[]{"TCATACTCGGAGGAGCTGG", "ACCTCCACCTCCCGG"}, null, false, 4, 4, 3});
        tests.add(new Object[]{"test3", null, null, new String[]{"TCATACTCGGAGGAGCTGG", "ACCTCCACCTCCCGG"}, false, 0, 0, null});

        tests.add(new Object[]{"test4", new String[]{"TCATACTCGGAGGAGCTGG", "ACCTCCACCTCCCGG"}, null, null, true, 4, 4, 3});
        tests.add(new Object[]{"test4", null, new String[]{"TCATACTCGGAGGAGCTGG", "ACCTCCACCTCCCGG"}, null, true, 4, 4, 3});
        tests.add(new Object[]{"test4", null, null, new String[]{"TCATACTCGGAGGAGCTGG", "ACCTCCACCTCCCGG"}, true, 0, 0, null});

        return tests.toArray(new Object[][]{});
    }

    private ArgumentsBuilder getBaseArgs(File fastqOut1, @Nullable File fastqOut2, @Nullable File summary) {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--fastq");
        args.add(normalizePath(getTestFile("fq1.fastq")));

        args.add("--output");
        args.add(normalizePath(fastqOut1));

        if (fastqOut2 != null) {
            args.add("--fastq2");
            args.add(normalizePath(getTestFile("fq2.fastq")));

            args.add("--output2");
            args.add(normalizePath(fastqOut2));
        }

        if (summary != null){
            args.add("--summaryFile");
            args.add(normalizePath(summary));
        }

        args.add("--tmp-dir");
        args.add(getTmpDir());

        return args;
    }

    @Test
    public void testEditDistanceBadInput() throws IOException {
        ArgumentsBuilder args = getBaseArgs(createTempFile("output", "-R1.fastq"), null, null);

        args.add("-e1");
        args.add("BADVALUES");

        args.add("--editDistance");
        args.add("1");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                0,
                UserException.BadInput.class);

        spec.executeTest("testEditDistance", this);

    }

    @Test
    public void testEditDistance() throws IOException {
        File output1 = createTempFile("output", "-R1.fastq");
        File output2 = createTempFile("output", "-R2.fastq");

        //Edit distance 1
        ArgumentsBuilder args = getBaseArgs(output1, output2, null);

        //One hit, ed=1
        args.add("-e1");
        args.add("TCTGCCTCTTG");

        args.add("-e");
        args.add("ACTGCTGCTTATT");

        args.add("--editDistance");
        args.add("1");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Collections.emptyList());

        spec.executeTest("testEditDistance1", this);

        IntegrationTestSpec.assertEqualTextFiles(getTestFile("ED1_R1.fastq"), output1);
        IntegrationTestSpec.assertEqualTextFiles(getTestFile("ED1_R2.fastq"), output2);

        //Edit distance 1
        args = getBaseArgs(output1, output2, null);

        args.add("-e1");
        args.add("TCTGCCTCTTG");

        args.add("-e");
        args.add("ACTGCTGCTTATT");

        args.add("--editDistance");
        args.add("2");

        spec = new IntegrationTestSpec(
                args.getString(),
                Collections.emptyList());

        spec.executeTest("testEditDistance2", this);

        IntegrationTestSpec.assertEqualTextFiles(getTestFile("ED2_R1.fastq"), output1);
        IntegrationTestSpec.assertEqualTextFiles(getTestFile("ED2_R2.fastq"), output2);
    }

    @Test(dataProvider = "testExpressionsData")
    public void testExpressionsPairedWithSummary(String testName, String[] exprs, String[] r1Exprs, String[] r2Exprs, boolean matchAllExpressions, int expectedLinesPE, int expectedLinesSE, @Nullable  Integer expectedSummaryLines) throws IOException {
        _testExpressionsPairedWithSummaryAndNames(testName, exprs, r1Exprs, r2Exprs, matchAllExpressions, true, false, expectedLinesPE, expectedSummaryLines);
    }

    @Test(dataProvider = "testExpressionsData")
    public void testExpressionsSingleEndWithSummary(String testName, String[] exprs, String[] r1Exprs, String[] r2Exprs, boolean matchAllExpressions, int expectedLinesPE, int expectedLinesSE, @Nullable  Integer expectedSummaryLines) throws IOException {
        _testExpressionsPairedWithSummaryAndNames(testName, exprs, r1Exprs, r2Exprs, matchAllExpressions, false, false, expectedLinesSE, expectedSummaryLines);
    }

    @Test(dataProvider = "testExpressionsData")
    public void testExpressionsPairedWithSummaryAndNames(String testName, String[] exprs, String[] r1Exprs, String[] r2Exprs, boolean matchAllExpressions, int expectedLinesPE, int expectedLinesSE, @Nullable  Integer expectedSummaryLines) throws IOException {
        _testExpressionsPairedWithSummaryAndNames(testName, exprs, r1Exprs, r2Exprs, matchAllExpressions, true, true, expectedLinesPE, expectedSummaryLines);
    }

    public void _testExpressionsPairedWithSummaryAndNames(String testName, String[] exprs, String[] r1Exprs, String[] r2Exprs, boolean matchAllExpressions, boolean paired, boolean addNames, int expectedLines, @Nullable  Integer expectedSummaryLines) throws IOException {
        File output1 = createTempFile("output", "-R1.fastq");

        File output2 = null;
        if (paired) {
            output2 = createTempFile("output", "-R2.fastq");
        }

        File summary = createTempFile("output", ".txt");

        ArgumentsBuilder args = getBaseArgs(output1, output2, summary);

        if (exprs != null) {
            Arrays.stream(exprs).forEach(x -> {
                args.add("-e");
                args.add(x);

                if (addNames) {
                    args.add("-en");
                    args.add("Name1");
                }
            });
        }

        if (r1Exprs != null) {
            Arrays.stream(r1Exprs).forEach(x -> {
                args.add("-e1");
                args.add(x);

                if (addNames) {
                    args.add("-e1n");
                    args.add("FName1");
                }
            });
        }

        if (paired && r2Exprs != null) {
            Arrays.stream(r2Exprs).forEach(x -> {
                args.add("-e2");
                args.add(x);

                if (addNames) {
                    args.add("-e2n");
                    args.add("RName1");
                }
            });
        }

        if (matchAllExpressions) {
            args.add("-ma");
        }

        runCommandLine(args);

        Assert.assertEquals(new XReadLines(output1).readLines().size(), expectedLines);
        if (output2 != null) {
            Assert.assertEquals(new XReadLines(output2).readLines().size(), expectedLines);
        }

        if (expectedSummaryLines == null) {
            expectedSummaryLines = 1 + (expectedLines / 4);
        }

        List<String> summaryLines = new XReadLines(summary).readLines();
        Assert.assertEquals(summaryLines.size(), expectedSummaryLines.intValue());

        if (addNames) {
            for (String line : summaryLines) {
                if (!line.startsWith("ReadName")) {
                    Assert.assertTrue(line.contains("Forward") ? line.contains("Name1") || line.contains("FName1") : line.contains("Name1") || line.contains("RName1"));
                }
            }
        }
        else {
            for (String line : summaryLines) {
                if (!line.startsWith("ReadName")) {
                    Assert.assertTrue(!line.contains("Name1"));
                }
            }
        }
    }
}
