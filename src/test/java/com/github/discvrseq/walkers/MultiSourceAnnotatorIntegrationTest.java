package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class MultiSourceAnnotatorIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/");

    @Test
    public void doBasicTest() throws Exception{
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-V");
        args.add(normalizePath(getInputVcf()));

        args.add("-cv");
        File clinvar = new File(testBaseDir, "walkers/MultiSourceAnnotator/clinvar.vcf");
        args.add(normalizePath(clinvar));
        ensureVcfIndex(clinvar);

        args.add("--liftoverReject");
        File liftover = new File(testBaseDir, "walkers/MultiSourceAnnotator/liftoverRejects.vcf");
        args.add(normalizePath(liftover));
        ensureVcfIndex(liftover);

        args.add("--cassandra");
        File cassandra = new File(testBaseDir, "walkers/MultiSourceAnnotator/cassandra.vcf");
        args.add(normalizePath(cassandra));
        ensureVcfIndex(cassandra);

        args.add("-O");
        args.add("%s");
        args.add("--tmp-dir");
        args.add(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(testBaseDir + "/walkers/MultiSourceAnnotator/multiSourceOutput.vcf"));

        spec.executeTest("doBasicTest", this);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "walkers/MultiSourceAnnotator/multiSourceInput.vcf");
    }
}
