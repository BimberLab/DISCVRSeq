package com.github.discvrseq.walkers.annotator;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;

public class MendelianViolationArgumentCollection {
    @Argument(fullName = "pedigreeValidationType", shortName = "pedValidationType", doc="The strictness for validating the pedigree.  Can be either STRICT or SILENT.  Default is STRICT", optional=true)
    public PedigreeValidationType pedigreeValidationType = PedigreeValidationType.STRICT;

    /**
     * This argument specifies the genotype quality (GQ) threshold that all members of a trio must have in order
     * for a site to be accepted as a mendelian violation. Note that the `-mv` flag must be set for this argument
     * to have an effect.
     */
    @Argument(fullName="mendelian-violation-qual-threshold", doc="Minimum GQ score for each trio member to accept a site as a violation", optional=true)
    public double minGenotypeQuality = 10.0;

    public static interface UsesMendelianViolationArgumentCollection {
        public void setArgumentCollection(MendelianViolationArgumentCollection args);
    }

    public void validateArguments() {

    }

    public SampleDB getSampleDB(GATKPath pedigreeFile) {
        return SampleDB.createSampleDBFromPedigree(pedigreeFile, pedigreeValidationType);
    }
}
