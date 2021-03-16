package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * This is a fairly specialized tool, designed to take two VCFs, and provide a summary table comparing their FILTER field by site. The input VCFs can be named, and these names will be used as
 * headers in the resulting TSV file.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar VcfFilterComparison \
 *     -R genomes.fasta \
 *     -V:WGS vcf1.vcf.gz \
 *     -V:WXS vcf2.vcf \
 *     -O output.txt
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This is a fairly specialized tool, designed to take two input VCFs, and provide a summary table comparing their FILTER field by site.",
        oneLineSummary = "Summarize filtering between two VCFs",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class VcfFilterComparison extends MultiVariantWalkerGroupedOnStart {
    @Argument(doc="File to which the summary should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        getDrivingVariantsFeatureInputs().forEach(fi -> drivingVariantNames.add(fi.getName()));
    }

    private List<String> drivingVariantNames = new ArrayList<>();

    private final Map<String, Long> combinations = new HashMap<>();

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        Map<String, String> vcfToFilters = new HashMap<>();

        for (FeatureInput<VariantContext> fi : getDrivingVariantsFeatureInputs()) {
            //This is very inefficient, but expedient until GATK incorporates a better mechanism to connect vc to source:
            List<VariantContext> vcs = features.getFeatures(fi, new SimpleInterval(variantContexts.get(0).getContig(), variantContexts.get(0).getStart(), variantContexts.get(0).getEnd()));

            List<String> filters = new ArrayList<>(vcs.stream().map(vc -> vc.isFiltered() ? vc.getFilters() : Collections.singleton("PASS")).flatMap(Collection::stream).collect(Collectors.toSet()));
            boolean isCalled = vcs.stream().map(vc -> vc.getNoCallCount() > 0).max(Boolean::compareTo).isPresent();

            vcfToFilters.put(fi.getName(), (isCalled ? "Y" : "N") + "|" + StringUtils.join(filters, ","));
        }

        String key = getKey(vcfToFilters);
        long val = combinations.getOrDefault(key, 0L);
        combinations.put(key, val + 1);
    }

    private String getKey(Map<String, String> vcfToFilters) {
        List<String> vals = new ArrayList<>();
        for (String name : drivingVariantNames)
        {
            vals.add(vcfToFilters.getOrDefault(name, "N|ND"));
        }

        return StringUtils.join(vals, "<>");
    }

    @Override
    public Object onTraversalSuccess() {
        try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(new File(outFile)), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            List<String> header = new ArrayList<>();
            drivingVariantNames.forEach(fi -> {
                header.add(fi);
                header.add(fi + "-Called");
            });
            writer.writeNext(header.toArray(new String[0]));

            for (String key : combinations.keySet()) {
                List<String> row = new ArrayList<>();
                for (String val : key.split("<>")) {
                    String[] vals = val.split(Pattern.quote("|"));
                    row.add(vals[1]);
                    row.add(vals[0]);
                }

                row.add(String.valueOf(combinations.get(key)));

                writer.writeNext(row.toArray(new String[0]));
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return null;
    }
}
