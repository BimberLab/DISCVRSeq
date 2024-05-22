package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;


/**
 * This walker will iterate a BAM and report sites where a high fraction of reads either begin or end soft clipping.
 * The original purpose is to identify positions that likely represent structural variants. This tool can be run with or without
 * a reference set of SVs.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar IdentifySoftClippedLoci \
 *     -R reference.fasta \
 *     -I bams.list \
 *     -O output.bed
 * </pre>
 * <p></p>
 *
 * Or you can provide a list of reference SVs. In this case, only sites that overlap the start/end of these SVs will be reported.
 * <pre>
 *  java -jar DISCVRseq.jar IdentifySoftClippedLoci \
 *     -I bams.list \
 *     -rv reference.vcf.gz \
 *     -O output.bed
 * </pre>
 * <p></p>
 * Would produce a BED file that looks like:
 * <p></p>
 * <pre>
 *     chr20 10000000 10000001 BeforeSoftClip   +   Sample1 0.5 100 200
 *     chr20 10000865 10000866 AfterSoftClip    +   Sample1 0.25 50 200    2-DEL-pbsv.DEL.10391.clr|2|132589349|INDEL
 * </pre>
 *
 * Where the columns are:
 * 1) Contig
 * 2) Start (0-based)
 * 3) End (0-based, exclusive)
 * 4) Reason for inclusion, either BeforeSoftClip or AfterSoftClip
 * 5) Strand (always listed as plus)
 * 6) Sample Name
 * 7) The percent of reads representing the start/end of the soft clip
 * 8) Total reads representing the start/end of the soft clip
 * 9) Total number of QC passing reads
 * 10) If a reference VCF is provided, and variants that overlap with this site, within --maxVariantDistance, will be reported
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool iterates a set of BAMs, reporting any loci (per sample) where the fraction of reads representing the start of end of a soft-clip is above the supplied threshold.",
        oneLineSummary = "Generate a list of loci with three or more distinct bases detected",
        programGroup = DiscvrSeqProgramGroup.class
)
public class IdentifySoftClippedLoci extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The output BED file", common = false, optional = true)
    private File outputFile = null;

    private PrintStream outputStream = null;

    @Advanced
    @Argument(fullName = "minDepth", shortName = "minDepth", doc = "Minimum QC+ read depth before a locus is considered callable", optional = true)
    int minDepth = 10;

    @Argument(fullName = "minFraction", shortName = "minFraction", doc = "Any site where at least this fraction of reads represent either the beginning or the end of a soft clip will be considered", optional = true)
    double minFraction = 0.1;

    @Argument(doc="Reference VCF", fullName = "refvcf", shortName = "rv", optional = true)
    public FeatureInput<VariantContext> referenceSites = null;

    @Argument(fullName = "maxVariantDistance", shortName = "mvd", doc = "If a reference VCF is provided, a loci is considered as overlapping the variant if the site is +/- this many bp from the start or end of the variant", optional = true)
    int variantDistance = 5;

    private enum REASON {
        BeforeSoftClip(),
        AfterSoftClip();
    }

    @Override
    public void onTraversalStart() {
        try {
            Utils.nonNull(outputFile);
            IOUtil.assertFileIsWritable(outputFile);
            outputStream = new PrintStream(outputFile);
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public void closeTool() {
        super.closeTool();

        outputStream.close();
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        assert alignmentContext.getStart() == alignmentContext.getEnd();

        List<VariantContext> vcs = new ArrayList<>();
        if (referenceSites != null) {
            vcs.addAll(featureContext.getValues(referenceSites));
        }

        for (SAMReadGroupRecord rg: getHeaderForReads().getReadGroups()) {
            ReadPileup pileup = alignmentContext.getBasePileup().getPileupForSample(rg.getSample(), getHeaderForReads());
            if (pileup.size() < minDepth) {
                continue;
            }

            AtomicInteger beforeSoftClip = new AtomicInteger();
            AtomicInteger afterSoftClip = new AtomicInteger();
            pileup.iterator().forEachRemaining(x -> {
                if (PileupCigarUtil.isAfterSoftOrHardClip(x)) {
                    afterSoftClip.getAndIncrement();
                }
                else if (PileupCigarUtil.isBeforeSoftOrHardClip(x)) {
                    beforeSoftClip.getAndIncrement();
                }
            });

            double beforePct = beforeSoftClip.get() / (double)pileup.size();
            if (beforePct > minFraction) {
                List<VariantContext> matching = vcs.stream().filter(vc -> Math.abs(vc.getStart() - alignmentContext.getStart()) < variantDistance).toList();
                printSite(alignmentContext, rg.getSample(), REASON.BeforeSoftClip, beforePct, pileup.size(), beforeSoftClip.get(), matching);
            }

            double afterPct = afterSoftClip.get() / (double)pileup.size();
            if (afterPct > minFraction) {
                List<VariantContext> matching = vcs.stream().filter(vc -> Math.abs(vc.getEnd() - alignmentContext.getStart()) < variantDistance).toList();
                printSite(alignmentContext, rg.getSample(), REASON.AfterSoftClip, afterPct, pileup.size(), afterSoftClip.get(), matching);
            }
        }
    }

    private void printSite(AlignmentContext alignmentContext, String sample, REASON reason, double pct, long depth, long totalReads, List<VariantContext> matching) {
        String vcs = matching.stream().map(vc -> {
            return new StringBuilder().append(vc.getID()).append("|").append(vc.getContig()).append("|").append(vc.getStart()).append("|").append(vc.getStructuralVariantType() == null ? vc.getType() : vc.getStructuralVariantType());
        }).collect(Collectors.joining(","));

        outputStream.printf("%s\t%d\t%d\t%s\t%s\t%s\t%f\t%d\t%d\t%s%n", alignmentContext.getContig(), alignmentContext.getStart()-1, alignmentContext.getEnd(), reason.name(), "+", sample, pct, depth, totalReads, vcs);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        List<ReadFilter> ret = super.getDefaultReadFilters();
        ret.add(new ReadFilterLibrary.NotDuplicateReadFilter());
        return ret;
    }

    // NOTE: this is directly copied from PileupElement. Consider trying to expose these methods directly in HTDJDK.
    public static class PileupCigarUtil
    {
        public static CigarOperator getAdjacentOperator(final PileupElement pileup, final Direction direction) {
            final int increment = direction.getIncrement();
            final int i = pileup.getCurrentCigarOffset() + increment;
            if (i < 0 || i >= pileup.getRead().numCigarElements()) {
                return null;
            }

            return pileup.getRead().getCigarElement(i).getOperator();
        }

        public static boolean isBeforeSoftOrHardClip(final PileupElement pileup)
        {
            return isImmediatelyBefore(pileup, CigarOperator.HARD_CLIP) || isImmediatelyBefore(pileup, CigarOperator.SOFT_CLIP);
        }

        public static boolean isAfterSoftOrHardClip(final PileupElement pileup)
        {
            return isImmediatelyAfter(pileup, CigarOperator.HARD_CLIP) || isImmediatelyAfter(pileup, CigarOperator.SOFT_CLIP);
        }

        private static boolean isImmediatelyAfter(final PileupElement pileup, final CigarOperator op) {
            return pileup.atStartOfCurrentCigar() && getAdjacentOperator(pileup, Direction.PREV) == op;
        }

        private static boolean isImmediatelyBefore(final PileupElement pileup, final CigarOperator op) {
            return pileup.atEndOfCurrentCigar() && getAdjacentOperator(pileup, Direction.NEXT) == op;
        }

        enum Direction { PREV(-1), NEXT(1);
            private final int increment;

            public int getIncrement() {
                return increment;
            }

            private Direction(final int increment){
                this.increment = increment;
            }
        }
    }
}