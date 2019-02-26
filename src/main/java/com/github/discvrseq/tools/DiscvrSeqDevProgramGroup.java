package com.github.discvrseq.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
 * Created by bimber on 8/1/2017.
 */
public class DiscvrSeqDevProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "Development/Internal Tools";
    }

    @Override
    public String getDescription() {
        return "Unlike the public tools group, these tools may have have received less rigorous testing and may be under active development.";
    }
}
