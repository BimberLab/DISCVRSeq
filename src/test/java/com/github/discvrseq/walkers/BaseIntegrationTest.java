package com.github.discvrseq.walkers;

import com.github.discvrseq.Main;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.BeforeClass;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.List;

public class BaseIntegrationTest extends CommandLineProgramTest {
    @Override
    @BeforeClass
    public void initGenomeLocParser() throws FileNotFoundException {
        //this expects files that dont exist
    }

    @Override
    public Object runCommandLine(final List<String> args) {
        //use our Main, not GATK's
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

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
        File expected = new File(input.getParent(), input.getName() + Tribble.STANDARD_INDEX_EXTENSION);
        if (expected.exists()){
            return;
        }

        Index index = IndexFactory.createDynamicIndex(input, new VCFCodec());
        try {
            index.writeBasedOnFeatureFile(input);
        }
        catch (IOException e){
            throw new RuntimeException(e);
        }
    }
}
