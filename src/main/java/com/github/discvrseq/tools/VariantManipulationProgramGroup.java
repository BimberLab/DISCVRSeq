package com.github.discvrseq.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

/**
 * Created by bimber on 8/1/2017.
 */
public class VariantManipulationProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "Variant Manipulation";
    }

    @Override
    public String getDescription() {
        return "These tools perform manipulation/QC of variant data, and should be rigorously validated and stable.";
    }
}
