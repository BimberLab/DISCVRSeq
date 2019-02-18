package com.github.discvrseq.walkers;

import com.github.discvrseq.Main;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.List;

public class BaseIntegrationTest extends BaseTest implements CommandLineProgramTester {
    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    public static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";

    private static final String publicTestDirRelative = "src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";

//    @Override
//    @BeforeClass
//    public void initGenomeLocParser() throws FileNotFoundException {
//        //this expects files that dont exist
//    }

    protected File downloadFile(String url, File saveTo) throws IOException {
        URL target = new URL(url);
        ReadableByteChannel rbc = Channels.newChannel(target.openStream());
        FileOutputStream fos = new FileOutputStream(saveTo);
        fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);

        return saveTo;
    }

    protected File downloadHg19Micro() throws IOException{
        File tmpDir = FileUtils.getTempDirectory();

        File fasta = downloadFile("https://raw.githubusercontent.com/broadinstitute/gatk/master/src/test/resources/hg19micro.fasta", new File(tmpDir, "hg19micro.fasta"));
        downloadFile("https://raw.githubusercontent.com/broadinstitute/gatk/master/src/test/resources/hg19micro.fasta.fai", new File(tmpDir, "hg19micro.fasta.fai"));
        downloadFile("https://raw.githubusercontent.com/broadinstitute/gatk/master/src/test/resources/hg19micro.dict", new File(tmpDir, "hg19micro.dict"));

        return fasta;
    }

    protected void ensureVcfIndex(File input){
        ensureIndex(input, new VCFCodec());
    }

    protected void ensureIndex(File input, FeatureCodec codec){
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
            tmpDir = normalizePath(IOUtils.getPath(System.getProperty("java.io.tmpdir")).toFile());
        }

        return tmpDir;
    }
}
