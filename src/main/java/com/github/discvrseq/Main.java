package com.github.discvrseq;

import com.google.common.annotations.VisibleForTesting;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by bimber on 7/30/2017.
 */
public final class Main extends org.broadinstitute.hellbender.Main {
    /** Only include discvr-seq tools . */
    @Override
    protected List<String> getPackageList() {
        final List<String> packageList = new ArrayList<>();
        packageList.add("com.github.discvrseq");
        return packageList;
    }

    /** Returns the command line that will appear in the usage. */
    protected String getCommandLineName() {
        return "java -jar DISCVRseq.jar";
    }

    public static void main(final String[] args) {
        final Main main = new Main();
        // if only the --version / -v is requested, exit directly
        if (main.printOnlyVersion(args)) {
            System.exit(0);
        }

        main.mainEntry(args);
    }

    @VisibleForTesting
    protected boolean printOnlyVersion(final String[] args) {
        System.out.println("Version: " + getClass().getPackage().getImplementationVersion());
        if (args.length == 1 && ("--version".equals(args[0]) || "-v".equals(args[0]))) {
            return true;
        }
        return false;
    }
}
