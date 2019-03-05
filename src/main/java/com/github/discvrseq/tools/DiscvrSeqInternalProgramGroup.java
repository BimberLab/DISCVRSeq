package com.github.discvrseq.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
 * Created by bimber on 8/1/2017.
 */
public class DiscvrSeqInternalProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "Specialized/Internal Tools";
    }

    @Override
    public String getDescription() {
        return "These tools are designed for fairly specialized purposes and are most likely only useful internally.";
    }
}
