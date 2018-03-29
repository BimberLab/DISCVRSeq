package com.github.discvrseq.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
 * Created by bimber on 8/1/2017.
 */
public class DiscvrSeqProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "DISCVR-Seq Tools (Public)";
    }

    @Override
    public String getDescription() {
        return "These tools should be rigorously validated and stable.";
    }
}
