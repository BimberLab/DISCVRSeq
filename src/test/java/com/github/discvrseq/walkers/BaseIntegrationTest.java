package com.github.discvrseq.walkers;

import com.github.discvrseq.Main;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class BaseIntegrationTest extends BaseTest implements CommandLineProgramTester {
    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    public static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";

    private static final String publicTestDirRelative = "src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";

    public static File testBaseDir = new File(publicTestDir + "com/github/discvrseq/TestData");

    protected File getHg19Micro() {
        return new File(testBaseDir, "hg19micro.fasta");
    }

    protected void ensureVcfIndex(File input){
        ensureIndex(input, new VCFCodec());
    }

    protected <FEATURE_TYPE extends Feature, SOURCE_TYPE> void ensureIndex(File input, FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec){
        File expected = new File(input.getParent(), input.getName() + Tribble.STANDARD_INDEX_EXTENSION);
        if (expected.exists()){
            return;
        }

        Index index = IndexFactory.createDynamicIndex(input, codec);
        try {
            index.writeBasedOnFeatureFile(input);
        }
        catch (IOException e){
            throw new RuntimeException(e);
        }
    }

    /**
     * This was added to so windows filepaths dont fail conversion to URIs in IOUtils
     * There must be a cleaner solution
     */
    public static String normalizePath(File file){
        return normalizePath(file.getPath());
    }

    public static String normalizePath(String path){
        path = path.replaceAll("\\\\+", "\\\\");
        int length = FilenameUtils.getPrefixLength(path);
        if (length > 0) {
            path = path.substring(length -1);
        }

        return FilenameUtils.separatorsToUnix(path);
    }

    /**
     * Returns the location of the resource directory. The default implementation points to the common directory for tools.
     */
    public static File getTestDataDir() {
        return new File("src/test/resources/com/github/discvrseq/");
    }

    @Override
    public String getTestedToolName() {
        return getTestedClassName();
    }

    @Override
    public Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

    @Override
    public Object runCommandLine(final List<String> args, final String toolName) {
        return new Main().instanceMain(makeCommandLineArgs(args, toolName));
    }

    private String tmpDir = null;

    public String getTmpDir() {
        if (tmpDir == null) {
            tmpDir = normalizePath(new File(System.getProperty("java.io.tmpdir")));
        }

        return tmpDir;
    }
}
