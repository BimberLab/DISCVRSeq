package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

public class MultiSourceAnnotatorIntegrationTest extends BaseIntegrationTest {
    private static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/");

    @Test
    public void doBasicTest() throws Exception{
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("-V");
        args.add(normalizePath(getInputVcf()));

        args.add("-O");
        File outFile = new File(normalizePath(getSafeNonExistentFile("multiSourceOutput.vcf")));
        args.add(outFile);

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

        runCommandLine(args.getArgsArray());

        File expected = new File(testBaseDir, "walkers/MultiSourceAnnotator/multiSourceOutput.vcf");
        IntegrationTestSpec.assertEqualTextFiles(new File(normalizePath(outFile)), expected);
    }

    private File getInputVcf(){
        return new File(testBaseDir, "walkers/MultiSourceAnnotator/multiSourceInput.vcf");
    }
}
