package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import com.github.discvrseq.util.CsvUtils;
import com.opencsv.ICSVWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * This produces a TSV report summarizing genotype qualities by Genotype Type.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar SummarizeGenotypeQuality \
 *     -V myFile.vcf \
 *     -O output.txt
 * </pre>
 *
 */
@CommandLineProgramProperties(summary="Tool for creating a summary of genotype qualities",
        oneLineSummary = "Tool for creating a summary of genotype qualities",
        programGroup = VariantManipulationProgramGroup.class)
public class SummarizeGenotypeQuality extends VariantWalker {
    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered sites", optional = true)
    public boolean excludeFiltered = true;

    private final Map<GenotypeType, Map<Integer, Long>> qualMap = new HashMap<>();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsWritable(new File(outFile));
    }

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (excludeFiltered && vc.isFiltered()) {
            return;
        }

        for (Genotype g : vc.getGenotypes()){
            if (excludeFiltered && g.isFiltered()) {
                continue;
            }

            GenotypeType t = g.getType();
            if (!qualMap.containsKey(t)) {
                qualMap.put(t, new HashMap<>());
            }

            int gq = g.getGQ();
            long count = qualMap.get(t).getOrDefault(gq, 0L);
            count++;
            qualMap.get(t).put(gq, count);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try (ICSVWriter writer = CsvUtils.getTsvWriter(new File(outFile))) {
            writer.writeNext(new String[]{"Type","Qual", "Total"});

            for (GenotypeType t : qualMap.keySet()) {
                for (int qual : qualMap.get(t).keySet()) {
                    writer.writeNext(new String[]{t.name(), String.valueOf(qual), String.valueOf(qualMap.get(t).get(qual))});
                }
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return super.onTraversalSuccess();
    }
}
