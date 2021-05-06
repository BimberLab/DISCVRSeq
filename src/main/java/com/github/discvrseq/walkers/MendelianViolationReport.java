package com.github.discvrseq.walkers;

import au.com.bytecode.opencsv.CSVWriter;
import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import com.github.discvrseq.walkers.annotator.MendelianViolationCount;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 *
 * This produces a TSV report summarizing the number of MVs detected per sample in the input VCF.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar MendelianViolationCount \
 *     -targetFasta targetGenomes.fasta \
 *     -ped pedigree.ped \
 *     -O output.txt
 * </pre>
 *
 */
@CommandLineProgramProperties(summary="Tool for creating a summary of MVs per sample from a VCF",
        oneLineSummary = "Tool for creating a summary of MVs per sample from a VCF",
        programGroup = DiscvrSeqProgramGroup.class)
public class MendelianViolationReport extends VariantWalker {
    @Argument(doc="File to which the report should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered sites", optional = true)
    private boolean excludeFiltered = true;

    @Argument(fullName="violationReportThreshold", shortName="rt", doc="Any sample with more than this many MVs will be reported", optional = true)
    private long violationReportThreshold = 500L;

    @Argument(fullName="contigsToSkip", shortName="cts", doc="A list of contigs to skip", optional = true)
    private List<String> contigsToSkip = new ArrayList<>();

    @Argument(fullName = StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName = StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for determining the population \"founders\"", optional=false)
    private GATKPath pedigreeFile;

    @Argument(fullName = "pedigreeValidationType", shortName = "pedValidationType", doc="The strictness for validating the pedigree.  Can be either STRICT or SILENT.  Default is STRICT", optional=true)
    private PedigreeValidationType pedigreeValidationType = PedigreeValidationType.STRICT;

    private Map<String, MVSummary> sampleMap;
    private SampleDB sampleDB = null;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        IOUtil.assertFileIsWritable(new File(outFile));

        sampleMap = new HashMap<>();
        VCFHeader header = getHeaderForVariants();
        List<String> samples = header.getGenotypeSamples();
        for (String sample : samples){
            sampleMap.put(sample, new MVSummary());
        }

        if ( sampleDB == null ) {
            sampleDB = initializeSampleDB();
        }
    }

    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(pedigreeValidationType);
        if (pedigreeFile != null)
            sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));

        if (!sampleMap.isEmpty()) {
            sampleDBBuilder.addSamplesFromSampleNames(sampleMap.keySet());
        }

        return sampleDBBuilder.getFinalSampleDB();
    }

    private static class MVSummary {
        long violationsDad = 0L;
        long violationsMom = 0L;
        long violationsTogether = 0L;
        long totalViolations = 0L;
        long totalCalled = 0L;

        public void addMV(MendelianViolationCount.MV mv, Genotype g){
            if (g.isCalled()) {
                totalCalled++;
            }

            if (mv == null){
                return;
            }

            if (mv.motherIsViolation){
                violationsMom++;
            }

            if (mv.fatherIsViolation){
                violationsDad++;
            }

            if (mv.violationCombined){
                violationsTogether++;
            }

            if (mv.isViolation()){
                totalViolations++;
            }
        }
    }

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (contigsToSkip != null && contigsToSkip.contains(referenceContext.getContig())) {
            return;
        }

        if (excludeFiltered && vc.isFiltered()) {
            return;
        }

        for (String sample : vc.getSampleNames()){
            Genotype g = vc.getGenotype(sample);
            if (g != null){
                Sample s = sampleDB.getSample(g.getSampleName());
                if (s != null){
                    MendelianViolationCount.MV mv = MendelianViolationCount.getMendelianViolation(s, vc, -1.0);
                    sampleMap.get(sample).addMV(mv, g);
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try (CSVWriter writer = new CSVWriter(IOUtil.openFileForBufferedUtf8Writing(new File(outFile)), '\t', CSVWriter.NO_QUOTE_CHARACTER)) {
            writer.writeNext(new String[]{"SampleName","TotalCalled","TotalViolations","MotherInconsistent","FatherInconsistent","InconsistentCombined","Mother","MotherHasData","MotherMVs","Father","FatherHasData","FatherMVs"});

            Set<String> samplesReported = new HashSet<>();
            Set<String> additionalSamplesToReport = new HashSet<>();
            for (String sn :sampleMap.keySet()) {
                MVSummary summary = sampleMap.get(sn);
                if (summary.totalViolations >= violationReportThreshold) {
                    samplesReported.add(sn);
                    reportSample(sn, summary, additionalSamplesToReport, writer);
                }
            }

            additionalSamplesToReport.removeAll(samplesReported);

            for(String sn :additionalSamplesToReport) {
                MVSummary summary = sampleMap.get(sn);
                if (summary != null) {
                    reportSample(sn, summary, Collections.emptySet(), writer);
                }
            }
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        return super.onTraversalSuccess();
    }

    private void reportSample(String sn, MVSummary summary, Set<String> additionalSamplesToReport, CSVWriter writer){
        Sample sample = sampleDB.getSample(sn);
        List<String> line = new ArrayList<>(Arrays.asList(
                sn,
                String.valueOf(summary.totalCalled),
                String.valueOf(summary.totalViolations),
                String.valueOf(summary.violationsMom),
                String.valueOf(summary.violationsDad),
                String.valueOf(summary.violationsTogether)
        ));

        appendParentToLine(sample.getMaternalID(), line, additionalSamplesToReport);
        appendParentToLine(sample.getPaternalID(), line, additionalSamplesToReport);

        writer.writeNext(line.toArray(new String[0]));
    }

    private void appendParentToLine(String parentId, List<String> line, Set<String> additionalSamplesToReport){
        if (parentId == null){
            line.add("Unknown");
            line.add("false");
            line.add("");
        }
        else {
            line.add(parentId);
            MVSummary summary = sampleMap.get(parentId);
            line.add(String.valueOf(summary != null));
            line.add(summary == null ? "" : String.valueOf(summary.totalViolations));
            additionalSamplesToReport.add(parentId);
        }
    }
}