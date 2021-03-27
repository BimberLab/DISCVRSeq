package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;

/**
 * Iterate a 10x BAM and summarize alignments per gene, including instances where reads map to multiple genes are are filtered.
 * The original purpose is to provide QC over 10x feature counting, such as identifying genes in a GTF that will be systematically under-counted.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar Summarize10xAlignments \
 *     -I possorted_bam.bam \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Iterate a 10x BAM and summarize alignments per gene",
        oneLineSummary = "Iterate a 10x BAM and summarize alignments per gene",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class Summarize10xAlignments extends ReadWalker {
    @Argument(doc="File to which the output table should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    Map<String, Long> readsPerGeneId = new HashMap<>();
    Map<String, Long> multimapReadsPerGeneId = new HashMap<>();
    Map<String, Set<String>> overlappingPerGeneId = new HashMap<>();

    Map<String, Long> readsPerGene = new HashMap<>();
    Map<String, Long> multimapReadsPerGene = new HashMap<>();
    Map<String, Set<String>> overlappingPerGene = new HashMap<>();

    SAMFileGATKReadWriter writer;

    @Override
    public void onTraversalStart() {
        IOUtil.assertFileIsWritable(new File(outFile));
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (read.hasAttribute("GX")) {
            appendValues(read, read.getAttributeAsString("GX"), readsPerGeneId, multimapReadsPerGeneId, overlappingPerGeneId);
        }

        if (read.hasAttribute("GN")) {
            appendValues(read, read.getAttributeAsString("GN"), readsPerGene, multimapReadsPerGene, overlappingPerGene);
        }
    }

    private void appendValues(GATKRead read, String attr, Map<String, Long> readsPerGene, Map<String, Long> multimapReadsPerGene, Map<String, Set<String>> overlappingPerGene) {
        Set<String> genes = new HashSet<>(Arrays.asList(attr.split(";")));
        genes.remove("None");
        if (genes.isEmpty()) {
            return;
        }

        if (genes.size() > 1) {
            for (String gene : genes) {
                multimapReadsPerGene.put(gene, multimapReadsPerGene.getOrDefault(gene, 0L) + 1);

                if (!overlappingPerGene.containsKey(gene)) {
                    overlappingPerGene.put(gene, new HashSet<>());
                }

                overlappingPerGene.get(gene).addAll(genes);
            }
        }
        else {
            readsPerGene.put(genes.iterator().next(), readsPerGene.getOrDefault(genes.iterator().next(), 0L) + 1);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        NumberFormat fmt = NumberFormat.getNumberInstance();
        fmt.setMaximumFractionDigits(2);

        try (CSVWriter csvWriter = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(new File(outFile)), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            csvWriter.writeNext(new String[]{"Type", "Gene", "TotalReads", "MultiMappedReads", "FractionMultiMapped", "OverlappingGenes"});

            Set<String> features = new TreeSet<>(readsPerGene.keySet());
            features.addAll(multimapReadsPerGene.keySet());
            for (String key : features) {
                Set<String> overlaps = overlappingPerGene.getOrDefault(key, Collections.emptySet());
                if (overlaps.contains(key)) {
                    overlaps.remove(key);
                }

                long denominator = readsPerGene.getOrDefault(key, 0L) + multimapReadsPerGene.getOrDefault(key, 0L);
                double pct = denominator == 0 ? 0 : (double)multimapReadsPerGene.getOrDefault(key, 0L) / denominator;
                csvWriter.writeNext(new String[]{"GeneName", key, String.valueOf(readsPerGene.getOrDefault(key, 0L)), String.valueOf(multimapReadsPerGene.getOrDefault(key, 0L)), fmt.format(pct), overlaps.isEmpty() ? "" : StringUtils.join(overlaps, ",")});
            }

            features = new TreeSet<>(readsPerGeneId.keySet());
            features.addAll(multimapReadsPerGeneId.keySet());
            for (String key : features) {
                Set<String> overlaps = overlappingPerGene.getOrDefault(key, Collections.emptySet());
                if (overlaps.contains(key)) {
                    overlaps.remove(key);
                }

                long denominator = readsPerGeneId.getOrDefault(key, 0L) + multimapReadsPerGeneId.getOrDefault(key, 0L);
                double pct = denominator == 0 ? 0 : (double)multimapReadsPerGeneId.getOrDefault(key, 0L) / denominator;
                csvWriter.writeNext(new String[]{"GeneId", key, String.valueOf(readsPerGeneId.getOrDefault(key, 0L)), String.valueOf(multimapReadsPerGeneId.getOrDefault(key, 0L)), fmt.format(pct), overlaps.isEmpty() ? "" : StringUtils.join(overlaps, ",")});
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return super.onTraversalSuccess();
    }
}
