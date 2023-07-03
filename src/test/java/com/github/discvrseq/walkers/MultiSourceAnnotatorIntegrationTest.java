package com.github.discvrseq.walkers;

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class MultiSourceAnnotatorIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/");

    @Test
    public void doBasicTest() throws Exception{
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("-V");
        args.addRaw(normalizePath(getInputVcf()));

        args.addRaw("-cv");
        File clinvar = new File(testBaseDir, "walkers/MultiSourceAnnotator/clinvar.vcf");
        args.addRaw(normalizePath(clinvar));
        ensureVcfIndex(clinvar);

        args.addRaw("--liftoverReject");
        File liftover = new File(testBaseDir, "walkers/MultiSourceAnnotator/liftoverRejects.vcf");
        args.addRaw(normalizePath(liftover));
        ensureVcfIndex(liftover);

        args.addRaw("--cassandra");
        File cassandra = new File(testBaseDir, "walkers/MultiSourceAnnotator/cassandra.vcf");
        args.addRaw(normalizePath(cassandra));
        ensureVcfIndex(cassandra);

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(testBaseDir + "/walkers/MultiSourceAnnotator/multiSourceOutput.vcf"));

        spec.executeTest("doBasicTest", this);
    }

    @Test
    public void doTestWithCustomFields() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addRaw("-V");
        args.addRaw(normalizePath(getInputVcf()));

        args.addRaw("-cv");
        File clinvar = new File(testBaseDir, "walkers/MultiSourceAnnotator/clinvar.vcf");
        args.addRaw(normalizePath(clinvar));
        ensureVcfIndex(clinvar);

        args.addRaw("--liftoverReject");
        File liftover = new File(testBaseDir, "walkers/MultiSourceAnnotator/liftoverRejects.vcf");
        args.addRaw(normalizePath(liftover));
        ensureVcfIndex(liftover);

        args.addRaw("--cassandra");
        File cassandra = new File(testBaseDir, "walkers/MultiSourceAnnotator/cassandra.vcf");
        args.addRaw(normalizePath(cassandra));
        ensureVcfIndex(cassandra);

        File cf = File.createTempFile("cassandra", ".args");
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(cf)) {
            writer.write("RFG\n");
            writer.write("FS_SN\n");
        }

        args.addRaw("-cf");
        args.addRaw(normalizePath(cf));

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(testBaseDir + "/walkers/MultiSourceAnnotator/multiSourceOutputSubsetArgs.vcf"));

        spec.executeTest("doTestWithCustomFields", this);

        cf.delete();
    }

    private File getInputVcf(){
        return new File(testBaseDir, "walkers/MultiSourceAnnotator/multiSourceInput.vcf");
    }
}
