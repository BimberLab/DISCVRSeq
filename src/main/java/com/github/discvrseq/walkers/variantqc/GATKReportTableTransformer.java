package com.github.discvrseq.walkers.variantqc;


import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.samples.SampleDB;

import javax.annotation.Nullable;

/**
 * Created by bimber on 5/31/2017.
 */
public interface GATKReportTableTransformer {
    public String getEvalModuleName();

    public GATKReportTable transform(GATKReportTable table, @Nullable SampleDB sampleDB);
}
