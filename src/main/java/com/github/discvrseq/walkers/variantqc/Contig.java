package com.github.discvrseq.walkers.variantqc;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.stream.Collectors;

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
            final SAMSequenceDictionary dict = getEngine().getSequenceDictionaryForDrivingVariants();
            if (dict == null) {
                throw new GATKException("Must provide a sequence dictionary, either for the genome FASTA or in the VCF header");
            }

            final Set<String> contigs = new LinkedHashSet<>();
            if (getEngine().getTraversalIntervals() == null) {
                dict.getSequences().stream().map(SAMSequenceRecord::getSequenceName).forEach(contigs::add);
            } else {
                getEngine().getTraversalIntervals().stream().map(SimpleInterval::getContig).forEach(contigs::add);
                contigs.forEach(contig -> {
                    if (dict.getSequence(contig) == null) {
                        throw new GATKException("The intervals used for this job contain a contig not present in the sequence dictionary: " + contig);
                    }
                });
            }

            if (contigs.size() > maxContigs) {
                // Sort based on sequence length:
                List<String> toRetain = contigs.stream().map(dict::getSequence).sorted(Collections.reverseOrder(Comparator.comparing(SAMSequenceRecord::getSequenceLength))).map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList());

                toRetain = toRetain.subList(0, maxContigs);

                if (contigsToRetain != null)
                {
                    for (String contig : contigsToRetain) {
                        if (!toRetain.contains(contig) && contigs.contains(contig)) {
                            toRetain.add(contig);
                        }
                    }
                }

                // Retain original order:
                List<String> finalSet = new ArrayList<>(contigs);
                finalSet.retainAll(toRetain);

                if (contigs.size() != finalSet.size()) {
                    finalSet.add(OTHER_CONTIG);
                }

                _contigNames = new LinkedHashSet<>(finalSet);
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
