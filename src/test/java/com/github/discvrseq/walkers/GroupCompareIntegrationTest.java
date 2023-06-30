package com.github.discvrseq.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class GroupCompareIntegrationTest extends BaseIntegrationTest {
    @Test
    public void doBasicTest() throws Exception {
        File fasta = getHg19Micro();
        ensureVcfIndex(getInputVcf());

        IntegrationTestSpec spec = new IntegrationTestSpec(
            " -R " + normalizePath(fasta) +
            " -V " + normalizePath(getInputVcf()) +
            " -G1 Sample1" +
            " -G1 Sample2" +
            " -O " + "%s" +
            " --tmp-dir " + getTmpDir(),
            Arrays.asList(getTestFile("groupCompareOutput.vcf").getPath()));

        spec.executeTest("doBasicTest", this);
    }

    @Test
    public void doBasicTestWithRef() throws Exception {
        File fasta = getHg19Micro();
        ensureVcfIndex(getInputVcf());
        ensureVcfIndex(getTestFile("groupCompareRef.vcf"));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -V " + normalizePath(getInputVcf()) +
                        " -RV " + normalizePath(getTestFile("groupCompareRef.vcf")) +
                        " -G1 Sample1" +
                        " -G1 Sample2" +
                        " -O " + "%s" +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(getTestFile("groupCompareOutputWithRef.vcf").getPath()));

        spec.executeTest("doBasicTest", this);
    }

    @Test
    public void doBasicTestWithRefAndSelect() throws Exception {
        File fasta = getHg19Micro();
        ensureVcfIndex(getInputVcf());
        ensureVcfIndex(getTestFile("groupCompareRef.vcf"));

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + normalizePath(fasta) +
                        " -V " + normalizePath(getInputVcf()) +
                        " -RV " + normalizePath(getTestFile("groupCompareRef.vcf")) +
                        " -G1 Sample1" +
                        " -G1 Sample2" +
                        " -O " + "%s" +
                        " -OT " + "%s" +
                        " -F REFFIELD" +
                        " -select " + "'G1_AF - REF_AF > 0.2'" +
                        " --tmp-dir " + getTmpDir(),
                Arrays.asList(
                        getTestFile("groupCompareOutputWithRef.vcf").getPath(),
                        getTestFile("groupCompareOutputWithRef.txt").getPath()
                ));

        spec.executeTest("doBasicTest", this);
    }

    private File getInputVcf() {
        return getTestFile("groupCompareInput1.vcf");
    }
}
