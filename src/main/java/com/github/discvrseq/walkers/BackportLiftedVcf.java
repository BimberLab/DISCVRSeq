package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqInternalProgramGroup;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This tool was originally created as part of an annotation pipeline for non-human data.  The input VCF from another species (or genome build) would be lifted to the human genome and annotated
 * in those coordinates. Lifeover must be performed using Picard Tools LiftoverVcf, which annotates lifted variants with the values for ORIGNAL_CONTIG, ORGINAL_START and ORIGINAL_ALLELE.  This tool
 * reads one of these lifted VCFs, and writes a new sorted VCF in which the original coordinates are restored.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar BackportLiftedVcf \
 *     -R currentGenome.fasta \
 *     -targetFasta targetGenomes.fasta \
 *     -V myVCF.vcf \
 *     -O output.vcf.gz
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This is a fairly specialized tool designed to restore the original coordinate in a VCF after liftover (created using Picard LiftoverVcf), back to the coordinates of the original genome.  It does this by reading the ORIGNAL_CONTIG, ORGINAL_START and ORIGINAL_ALLELE annotations left by Picard.",
        oneLineSummary = "Backport lifted VCF to the original coordinates",
        programGroup = DiscvrSeqInternalProgramGroup.class
)
public class BackportLiftedVcf extends VariantWalker {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String outFile = null;

    @Argument(doc="Target FASTA file", fullName = "targetFasta", shortName = "targetFasta", optional = false)
    public String targetGenome = null;

    @Argument(doc = "The liftOver can result in dropped ALT alleles. If this flag is provided, the tool will error if a given loci does not have the same number of ref/alt alleles. The default is to skip it.",  fullName = "errorOnLociWithMismatchedAlleles", shortName = "errorOnLociWithMismatchedAlleles", optional = true)
    public boolean errorOnLociWithMismatchedAlleles = false;

    private SortingCollection<VariantContext> sorter;
    private VCFHeader outputHeader;

    private final int MAX_RECORDS_IN_RAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();

    public static final String ORIGINAL_CONTIG = "OriginalContig";
    public static final String ORIGINAL_START = "OriginalStart";
    public static final String ORIGINAL_ALLELES = "OriginalAlleles";

    public static final String ORIGINAL_CONTIG_BCF = "SRC_CHROM";
    public static final String ORIGINAL_START_BCF = "SRC_POS";
    public static final String ORIGINAL_ALLELES_BCF = "SRC_REF_ALT";

    public static final String LIFTED_CONTIG = "LiftedContig";
    public static final String LIFTED_START = "LiftedStart";
    public static final String LIFTED_STOP = "LiftedStop";

    private static final List<VCFInfoHeaderLine> ATTRS = CollectionUtil.makeList(
            new VCFInfoHeaderLine(LIFTED_CONTIG, 1, VCFHeaderLineType.String, "The name of the contig/chromosome after liftover."),
            new VCFInfoHeaderLine(LIFTED_START, 1, VCFHeaderLineType.Integer, "The start position of the variant on the lifted contig."),
            new VCFInfoHeaderLine(LIFTED_STOP, 1, VCFHeaderLineType.Integer, "The end position of the variant on the lifted contig.")
    );

    private long skippedLoci = 0;

    private final Log log = Log.getInstance(BackportLiftedVcf.class);

    @Override
    public void onTraversalStart() {
        Utils.nonNull(outFile);

        if (getHeaderForVariants().getInfoHeaderLine(ORIGINAL_CONTIG) == null && getHeaderForVariants().getInfoHeaderLine(ORIGINAL_CONTIG_BCF) == null)
        {
            throw new IllegalArgumentException("The input VCF lacks both the " + ORIGINAL_CONTIG + " and " + ORIGINAL_CONTIG_BCF + " annotations");
        }

        if (getHeaderForVariants().getInfoHeaderLine(ORIGINAL_START) == null && getHeaderForVariants().getInfoHeaderLine(ORIGINAL_START_BCF) == null)
        {
            throw new IllegalArgumentException("The input VCF lacks both the " + ORIGINAL_START + " and " + ORIGINAL_START_BCF + " annotations");
        }

        if (getHeaderForVariants().getInfoHeaderLine(ORIGINAL_ALLELES) == null && getHeaderForVariants().getInfoHeaderLine(ORIGINAL_ALLELES_BCF) == null)
        {
            throw new IllegalArgumentException("The input VCF lacks both the " + ORIGINAL_ALLELES + " and " + ORIGINAL_ALLELES_BCF + " annotations");
        }

        prepareVcfHeader();
        initializeSorter(outputHeader);
    }

    private void prepareVcfHeader(){
        VCFHeader header = getHeaderForVariants();

        Set<VCFHeaderLine> lines = new LinkedHashSet<>(header.getMetaDataInInputOrder());
        lines.remove(header.getInfoHeaderLine(ORIGINAL_CONTIG));
        lines.remove(header.getInfoHeaderLine(ORIGINAL_START));
        lines.remove(header.getInfoHeaderLine(ORIGINAL_ALLELES));
        lines.remove(header.getInfoHeaderLine(ORIGINAL_CONTIG_BCF));
        lines.remove(header.getInfoHeaderLine(ORIGINAL_START_BCF));
        lines.remove(header.getInfoHeaderLine(ORIGINAL_ALLELES_BCF));

        outputHeader = new VCFHeader(lines, header.getSampleNamesInOrder());

        File targetFasta = new File(targetGenome);
        IOUtil.assertFileIsReadable(targetFasta);

        File dict = new File(ReferenceUtils.getFastaDictionaryFileName(targetGenome));
        IOUtil.assertFileIsReadable(dict);

        SAMSequenceDictionary sequenceDictionary = ReferenceUtils.loadFastaDictionary(dict);
        outputHeader.setSequenceDictionary(sequenceDictionary);

        for (final VCFInfoHeaderLine line : ATTRS) {
            outputHeader.addMetaDataLine(line);
        }
    }

    private void initializeSorter(VCFHeader outputHeader) {
        File tmpDir = IOUtil.getDefaultTmpDir();
        if (!tmpDir.exists()) {
            tmpDir.mkdirs();
        }

        sorter = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(outputHeader, true),
                outputHeader.getVCFRecordComparator(),
                MAX_RECORDS_IN_RAM, tmpDir.toPath());
    }

    private String getCoalescedFieldName(VariantContext variant, String fn1, String fn2)
    {
        if (variant.hasAttribute(fn1))
        {
            return fn1;
        }
        else if (variant.hasAttribute(fn2))
        {
            return fn2;
        }

        throw new IllegalStateException("Variant lacks " + fn1 + " or " + fn2 + " fields. This should have been caught during initialization");
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        String origChr = variant.getContig();
        int origStart = variant.getStart();
        int origStop = variant.getEnd();

        String targetChr = variant.getAttributeAsString(getCoalescedFieldName(variant, ORIGINAL_CONTIG, ORIGINAL_CONTIG_BCF), null);
        Utils.nonNull(targetChr, "Missing annotation for " + ORIGINAL_CONTIG + " at position: " + origChr + " " + origStart);

        if (outputHeader.getSequenceDictionary().getSequence(targetChr) == null){
            throw new IllegalArgumentException("Unknown contig: " + targetChr + " at position: " + origChr + " " + origStart);
        }

        int targetStart = variant.getAttributeAsInt(getCoalescedFieldName(variant, ORIGINAL_START, ORIGINAL_START_BCF), -1);
        if (targetStart == -1){
            throw new IllegalArgumentException("Missing annotation for " + ORIGINAL_START + " at position: " + origChr + " " + origStart);
        }

        List<Allele> targetAlleles = new ArrayList<>();
        List<String> origAlleles = variant.getAttributeAsStringList(getCoalescedFieldName(variant, ORIGINAL_ALLELES, ORIGINAL_ALLELES_BCF), null);
        if (origAlleles == null || origAlleles.isEmpty()){
            targetAlleles = variant.getAlleles();
        }
        else {
            if (origAlleles.size() != variant.getAlleles().size()) {
                if (errorOnLociWithMismatchedAlleles) {
                    throw new IllegalArgumentException("Original alleles listed for position: " + origChr + " " + origStart + " do not have the same number as the alleles at this site. Existing: " + variant.getAlleles().stream().map(Allele::getDisplayString).collect(Collectors.joining(",")) + ", annotation: " + origAlleles.stream().collect(Collectors.joining(",")));
                }

                skippedLoci++;
                return;
            }

            int idx = 0;
            for (String bases : origAlleles){
                targetAlleles.add(Allele.create(bases, idx ==0));
                idx++;
            }
        }

        //the original variant is always in forward orientation
        int targetEnd = targetStart + targetAlleles.get(0).length() - 1;

        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        vcb.chr(targetChr);
        vcb.start(targetStart);
        vcb.stop(targetEnd);
        vcb.alleles(targetAlleles);

        vcb.genotypes(fixGenotypes(variant.getGenotypes(), variant.getAlleles(), vcb.getAlleles()));

        vcb.attribute(LIFTED_CONTIG, origChr);
        vcb.attribute(LIFTED_START, origStart);
        vcb.attribute(LIFTED_STOP, origStop);

        vcb.rmAttribute(ORIGINAL_CONTIG);
        vcb.rmAttribute(ORIGINAL_START);
        vcb.rmAttribute(ORIGINAL_ALLELES);

        vcb.rmAttribute(ORIGINAL_CONTIG_BCF);
        vcb.rmAttribute(ORIGINAL_START_BCF);
        vcb.rmAttribute(ORIGINAL_ALLELES_BCF);

        sorter.add(vcb.make());
    }

    private GenotypesContext fixGenotypes(final GenotypesContext originals, List<Allele> originalAlleles, List<Allele> newAlleles) {
        if (originalAlleles.equals(newAlleles)) {
            return originals;
        }

        final GenotypesContext fixedGenotypes = GenotypesContext.create(originals.size());
        for (final Genotype genotype : originals) {
            final List<Allele> fixedAlleles = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                if (allele.isSymbolic() || allele.isNoCall()){
                    fixedAlleles.add(allele);
                }
                else {
                    int idx = originalAlleles.indexOf(allele);
                    if (idx == -1) {
                        throw new IllegalStateException("Allele not found: " + allele.toString() + ", " + originalAlleles + "/ " + newAlleles);
                    }
                    fixedAlleles.add(newAlleles.get(idx));
                }
            }
            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
        }
        return fixedGenotypes;
    }

    @Override
    public Object onTraversalSuccess() {
        writeSortedOutput();

        log.info("Total loci skipped due to mismatched allele number: " + skippedLoci);

        return null;
    }

    private void writeSortedOutput() {
        final ProgressLogger writeProgress = new ProgressLogger(log, 25000, "wrote", "records");
        final EnumSet<Options> options = createOutputVariantIndex ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
        if (lenientVCFProcessing) {
            options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        }

        try (final VariantContextWriter out = GATKVariantContextUtils.createVCFWriter(
                new File(outFile).toPath(),
                outputHeader.getSequenceDictionary(),
                createOutputVariantMD5,
                options.toArray(new Options[options.size()]));CloseableIterator<VariantContext> it = sorter.iterator()) {

            out.writeHeader(outputHeader);
            while (it.hasNext())
            {
                VariantContext variantContext = it.next();
                out.add(variantContext);
                writeProgress.record(variantContext.getContig(), variantContext.getStart());
            }
        }
        finally {
            sorter.cleanup();
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }
}
