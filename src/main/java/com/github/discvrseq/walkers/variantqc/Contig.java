package com.github.discvrseq.walkers.variantqc;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;

/**
 * Stratifies the evaluation by each contig in the reference sequence. Note: if the user supplies custom intervals, it will defer to these rather than the full sequence dictionary
 */
public class Contig extends VariantStratifier {
    private final String OTHER_CONTIG = "Other";

    private Set<String> _contigNames = null;

    public Contig(VariantEvalEngine engine, int maxContigs, List<String> contigsToRetain) {
        super(engine);

        states.addAll(getContigNames(maxContigs, contigsToRetain));
        states.add("all");
    }

    /**
     * @return The list of contig names to be traversed, preferentially taking user supplied intervals, but otherwise defaulting to driving variants
     */
    private Set<String> getContigNames(int maxContigs, List<String> contigsToRetain) {
        if (_contigNames == null) {
            final Set<String> contigs = new LinkedHashSet<>();
            if (getEngine().getTraversalIntervals() == null) {
                getEngine().getSequenceDictionaryForDrivingVariants().getSequences().stream().sorted(Collections.reverseOrder(Comparator.comparing(SAMSequenceRecord::getSequenceLength))).map(SAMSequenceRecord::getSequenceName).forEach(contigs::add);
            } else {
                getEngine().getTraversalIntervals().stream().map(SimpleInterval::getContig).forEach(contigs::add);
            }

            if (contigs.size() > maxContigs) {
                List<String> subset = new ArrayList<>(contigs);
                subset = subset.subList(0, maxContigs);

                if (contigsToRetain != null)
                {
                    for (String contig : contigsToRetain) {
                        if (!subset.contains(contig) && contigs.contains(contig)) {
                            subset.add(contig);
                        }
                    }
                }

                if (contigs.size() != subset.size()) {
                    subset.add(OTHER_CONTIG);
                }

                _contigNames = new LinkedHashSet<>(subset);
            }
            else {
                _contigNames = contigs;
            }
        }

        return _contigNames;
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null) {
            return Arrays.asList("all", _contigNames.contains(eval.getContig()) ? eval.getContig() : OTHER_CONTIG);
        } else {
            return Collections.emptyList();
        }
    }
}
