package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;


/**
 * This tool will compare the samples in a VCF against and number of reference VCFs.  It will produce a report summarizing the number of alleles shared between each sample and the reference sets
 * One possible use of this tool is to compare VCF samples against reference panels of alleles, such as alleles unique to different ethnicities or populations.
 *
 * The format of the reference VCF(s) needs to be somewhat specific.  Each site must be either:
 *  - A single allele, in which case it will score any sample containing that allele
 *  - Two alleles, in which case it will score any sample containing the ALT allele
 *
 *  Filtered sites on the reference VCFs are not allowed.  Also be aware that this tools assumes the size of the reference files is
 *  relatively small (though in theory 100Ks might work), since the intervals are read into memory.
 *
 * <h3>Usage example (note naming of the --ref-sites files):</h3>
 * <pre>
 *  java -jar DISCVRseq.jar VariantConcordanceScore \
 *     --ref-sites:SET1 ref1.vcf \
 *     --ref-sites:SET2 ref2.vcf.gz \
 *     -V myVCF.vcf.gz \
 *     -O output.txt
 * </pre>
 *
 * This will output a table with one row for each sample/reference pair, with the following data:
 * <ul>
 *  <li>SampleName: name of sample</li>
 *  <li>ReferenceName: name of the reference, as given on the command line ("i.e. --ref-sites:SET1 ref1.vcf)</li>
 *  <li>MarkersMatched: The number of sites from this sample where an allele matched the allele from this reference</li>
 *  <li>MarkersMismatched: The number of sites from this sample where there was callable data but no alleles matched the allele from this reference</li>
 *  <li>FractionMatched: The fraction of markers that matched (which does not count without callable data in this sample)</li>
 *  <li>TotalMarkersForSet: The total markers in this reference set</li>
 *  <li>FractionWithData: The fraction of markers from this reference set that had callable genotypes in this sample</li>
 * </ul>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool compares an input VCF against a set of reference VCFs, reporting concordance between them",
        oneLineSummary = "Summarize allele-level concordance between a VCF and reference VCFs",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class VariantConcordanceScore extends VariantWalker {

    @Argument(fullName = "ref-sites", shortName = "rs", doc = "VCF file containing sites to test.  Must be uniquely named", optional = false)
    public List<FeatureInput<VariantContext>> referenceFiles = new ArrayList<>();

    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    private Map<SimpleInterval, Map<Allele, Set<String>>> refMap = null;
    private Map<String, Long> totalMarkerByRef = new HashMap<>();

    @Override
    public List<SimpleInterval> getTraversalIntervals() {
        if (refMap == null) {
            prepareRefIntervals();
        }

        List<SimpleInterval> ret = new ArrayList<>(refMap.keySet());
        Collections.sort(ret, new Comparator<SimpleInterval>() {
            @Override
            public int compare(SimpleInterval o1, SimpleInterval o2) {
                if (o2 == null) return -1; // nulls last

                int result = o1.getContig().compareTo(o2.getContig());
                if (result == 0) {
                    if (o1.getStart() == o2.getStart()) {
                        result = o1.getEnd() - o2.getEnd();
                    } else {
                        result = o1.getStart() - o2.getStart();
                    }
                }

                return result;
            }
        });

        return ret;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsWritable(new File(outFile));

        prepareRefIntervals();
    }

    private void prepareRefIntervals() {
        Map<SimpleInterval, Map<Allele, Set<String>>> ret = new HashMap<>();

        referenceFiles.stream().forEach(
            f -> {
                try (VCFFileReader reader = new VCFFileReader(f.toPath()); CloseableIterator<VariantContext> it = reader.iterator()) {
                    while (it.hasNext()) {
                        VariantContext vc = it.next();
                        if (vc.isFiltered()) {
                            throw new IllegalStateException("Reference VCFs should not have filtered variants: " + vc.toStringWithoutGenotypes());
                        }

                        if (vc.getAlternateAlleles().size() > 1) {
                            throw new IllegalStateException("Reference has site with multiple alternates.  Must have either zero or one ALT: " + vc.toStringWithoutGenotypes());
                        }

                        SimpleInterval i = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());
                        Allele a = vc.isVariant() ? vc.getAlternateAllele(0) : vc.getReference();
                        Map<Allele, Set<String>> map = ret.getOrDefault(i, new HashMap<>());
                        Set<String> sources = map.getOrDefault(a, new HashSet<>());

                        //NOTE: should we throw warning or error if this is >1, indicating references are not mutually exclusive?
                        sources.add(f.getName());
                        map.put(a, sources);
                        ret.put(i, map);

                        totalMarkerByRef.put(f.getName(), totalMarkerByRef.getOrDefault(f.getName(), 0L) + 1);
                    }
                }
            }
        );

        refMap = ret;
        logger.info("total intervals: " + refMap.size());
    }

    private Map<String, SampleStats> sampleMap = new HashMap<>();

    private class SampleStats {
        long totalNoCall = 0;
        Map<String, Long> hits = new HashMap<>();
        Map<String, Long> misses = new HashMap<>();
    }

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (vc.isFiltered()) {
            return;
        }

        SimpleInterval i = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());
        if (refMap.containsKey(i)) {
            Map<Allele, Set<String>> map = refMap.get(i);

            for (Genotype g : vc.getGenotypes()) {
                SampleStats ss = sampleMap.containsKey(g.getSampleName()) ? sampleMap.get(g.getSampleName()) : new SampleStats();
                if (!g.isCalled()) {
                    ss.totalNoCall += 1;
                    continue;
                }

                for (Allele a : map.keySet()) {
                    if (g.getAlleles().contains(a)) {
                        for (String refHit : map.get(a)) {
                            long total = ss.hits.getOrDefault(refHit, 0L);
                            total++;
                            ss.hits.put(refHit, total);
                        }
                    }
                    else {
                        for (String refHit : map.get(a)) {
                            long total = ss.misses.getOrDefault(refHit, 0L);
                            total++;
                            ss.misses.put(refHit, total);
                        }
                    }
                }

                sampleMap.put(g.getSampleName(), ss);
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        //Note: need should consider implementing some kind of logic to make actual calls per reference
        NumberFormat format = NumberFormat.getInstance();
        format.setMaximumFractionDigits(2);

        try (CSVWriter output = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(new File(outFile)), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            output.writeNext(new String[]{"SampleName", "ReferenceName", "MarkersMatched", "MarkersMismatched", "FractionMatched", "TotalMarkersForSet", "FractionWithData"});

            for (String sample : sampleMap.keySet()) {
                SampleStats ss = sampleMap.get(sample);
                for (String ref : totalMarkerByRef.keySet()) {
                    long totalWithData = ss.hits.getOrDefault(ref, 0L) + ss.misses.getOrDefault(ref, 0L);
                    double fractionCalled = totalWithData == 0 ? 0 : (double)ss.hits.getOrDefault(ref, 0L) / totalWithData;
                    double fractionWithData = totalMarkerByRef.get(ref) == 0 ? 0 : (double)totalWithData /  totalMarkerByRef.get(ref);

                    output.writeNext(new String[]{sample, ref, String.valueOf(ss.hits.getOrDefault(ref, 0L)), String.valueOf(ss.misses.getOrDefault(ref, 0L)), format.format(fractionCalled), String.valueOf(totalMarkerByRef.get(ref)), format.format(fractionWithData)});
                }
            }

        }
        catch (IOException e) {
            throw new GATKException(e.getMessage());
        }

        return super.getTraversalIntervals();
    }
}
