package com.github.discvrseq;

import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;

public class TestUtils {
    // this was added to so windows filepaths dont fail conversion to URIs in IOUtils
    // there must be a cleaner solution
    public static String fixFilePath(File file){
        file = IOUtils.absolute(file);
        String ret = file.toURI().toString();
        ret = ret.replaceAll("file://", "");
        ret = ret.replaceAll("//", "/");

        return ret;
    }
}
