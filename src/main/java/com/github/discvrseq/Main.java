package com.github.discvrseq;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Log;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by bimber on 7/30/2017.
 */
public final class Main extends org.broadinstitute.hellbender.Main {
    /** Stream output for exceptions. */
    @VisibleForTesting
    PrintStream exceptionOutput = System.err;

    /** Stream output for results. */
    @VisibleForTesting
    PrintStream resultOutput = System.out;

    /** Only include discvr-seq tools . */
    @Override
    protected List<String> getPackageList() {
        final List<String> packageList = new ArrayList<>();
        packageList.add("com.github.discvrseq");
        return packageList;
    }

    /** Note: currently no single class is included. */
    @Override
    protected List<Class<? extends CommandLineProgram>> getClassList() {
        // TODO: explore other tools from the GATK4 framework that may be useful for our toolkit
        return Collections.emptyList();
    }

    /** Returns the command line that will appear in the usage. */
    protected String getCommandLineName() {
        return "java -jar DiscvrAnalysisToolkit.jar";
    }

    public static void main(final String[] args) {
        final Main main = new Main();
        // if only the --version / -v is requested, exit directly
        if (main.printOnlyVersion(args)) {
            System.exit(0);
        }
        new Main().mainEntry(args);
    }

    @VisibleForTesting
    protected boolean printOnlyVersion(final String[] args) {
        if (args.length == 1 && ("--version".equals(args[0]) || "-v".equals(args[0]))) {
            handleResult(getClass().getPackage().getImplementationVersion());
            return true;
        }
        return false;
    }

    @Override
    protected void handleResult(final Object result) {
        // TODO: print something else and/or handle metrics?
        if (result != null) {
            resultOutput.println(result);
        }
    }

    @Override
    protected void handleUserException(Exception e) {
        printDecoratedExceptionMessage(exceptionOutput, e, "A USER ERROR has occurred: ");
        printStackTrace(e);
    }

    @Override
    protected void handleNonUserException(final Exception e) {
        printDecoratedExceptionMessage(exceptionOutput, e, "UNEXPECTED ERROR: ");
        exceptionOutput
                .println("Please, search for this error in our issue tracker or post a new one:");
        printStackTrace(e);
    }

    /**
     * Prints the stack-trace into {@link #exceptionOutput} only if
     * {@link htsjdk.samtools.util.Log.LogLevel#DEBUG} is enabled.
     */
    @VisibleForTesting
    protected final void printStackTrace(final Exception e) {
        if (Log.isEnabled(Log.LogLevel.DEBUG)) {
            e.printStackTrace(exceptionOutput);
        }
    }
}
