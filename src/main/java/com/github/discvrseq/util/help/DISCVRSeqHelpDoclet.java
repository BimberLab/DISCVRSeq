package com.github.discvrseq.util.help;

import jdk.javadoc.doclet.DocletEnvironment;
import org.broadinstitute.hellbender.utils.help.GATKHelpDoclet;

@SuppressWarnings("removal")
public class DISCVRSeqHelpDoclet extends GATKHelpDoclet {
    public DISCVRSeqHelpDoclet() {

    }

    public static boolean processDocs(final DocletEnvironment docletEnv) {
        return new DISCVRSeqHelpDoclet().run(docletEnv);
    }
}
