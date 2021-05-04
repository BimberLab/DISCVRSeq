package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.DiscvrSeqProgramGroup;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * MergeVcfsAndGenotypes is an adaptation of the Broad Institute GATK 3 tool CombineVariants. This is different from Picard
 * tools MergeVcfs in that it will not only merge sites, but also merge genotypes.
 *
 * <ul>
 * <li><b>MERGE:</b> combines multiple variant records present at the same site in the different input sources into a
 * single variant record in the output. If sample names overlap, then they are "uniquified" by default, which means a
 * suffix is appended to make them unique. <em>Note that in version 3.3, the automatic uniquifying was disabled
 * (unintentionally), and required setting `-genotypeMergeOption UNIQUIFY` manually.</em></li>
 *
 * <li><b>UNION:</b> assumes that each ROD source represents the same set of samples (although this is not enforced).
 * It uses the priority list (if provided) to emit a single record instance at every position represented in the input RODs.</li>
 * </ul>
 *
 * <p>By default, the input sets will be named variants, variants2, variants3, and so on. You can override this by
 * providing an explicit name tag for each input, using the syntax " -V:name,format". Each input tagged in this
 * way will be labeled as such in the output (i.e., set=name rather than set=variants2). For example, you could specify
 * a set of control samples as " -V:control,vcf my_control_samples.vcf", and the resulting VCF records would contain
 * the annotation "set=control" in the INFO field. It is strongly recommended to provide explicit names in this way
 * when a rod priority list is provided.</p>
 *
 * <p>MergeVcfsAndGenotypes will emit a record for every site that was present in any of your input VCF files, and will annotate
 * (in the set attribute in the INFO field) whether the record had a PASS or FILTER status in each input ROD . In effect,
 * MergeVcfsAndGenotypes always produces a union of the input VCFs.  However, any part of the Venn of the merged VCFs
 * can be extracted using JEXL expressions on the set attribute using SelectVariants.  If you want to extract just
 * the records in common between two VCFs, you would first run MergeVcfsAndGenotypes on the two files to generate a single
 * VCF and then run SelectVariants to extract the common records with `-select 'set == "Intersection"'`, as worked out
 * in the detailed example in the documentation guide.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more variant sets to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF.
 * </p>
 *
 * <h3>Usage examples</h3>
 * &nbsp;
 * <h4>Merge two separate callsets</h4>
 * <pre>
 * java -jar DISCVRseq.jar \
 *   -T MergeVcfsAndGenotypes \
 *   -R reference.fasta \
 *   --variant input1.vcf \
 *   --variant input2.vcf \
 *   -O output.vcf \
 *   -genotypeMergeOption UNIQUIFY
 * </pre>
 *
 * <h4>Get the union of calls made on the same samples </h4>
 * <pre>
 * java -jar DISCVRseq.jar \
 *   -T MergeVcfsAndGenotypes \
 *   -R reference.fasta \
 *   --variant:foo input1.vcf \
 *   --variant:bar input2.vcf \
 *   -O output.vcf \
 *   -genotypeMergeOption PRIORITIZE \
 *   -priority foo,bar
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>This tool is not intended to manipulate gVCFS! To combine GVCF files output for different samples by HaplotypeCaller, use CombineGVCFs.</li>
 * </ul>
 *
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This will merge two or more input VCF files at the site and genotype level.",
        oneLineSummary = "Generate a merged VCF from one or more inputs",
        programGroup = DiscvrSeqProgramGroup.class
)
public class MergeVcfsAndGenotypes extends MultiVariantWalkerGroupedOnStart {
    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public File out = null;

    @Argument(shortName="genotypeMergeOption", doc="Determines how we should merge genotype records for samples shared across the VCF files", optional = true)
    public GATKVariantContextUtils.GenotypeMergeType genotypeMergeOption = null;

    @Argument(shortName="filteredRecordsMergeType", doc="Determines how we should handle records seen at the same site in the VCF, but with different FILTER fields", optional = true)
    public GATKVariantContextUtils.FilteredRecordMergeType filteredRecordsMergeType = GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED;

    @Hidden
    @Argument(shortName="multipleAllelesMergeType", doc="Determines how we should handle records seen at the same site in the VCF, but with different allele types (for example, SNP vs. indel)", optional = true)
    public MultipleAllelesMergeType multipleAllelesMergeType = MultipleAllelesMergeType.BY_TYPE;

    /**
     * Refers to the merging priority behavior described in the tool documentation regarding the choice of which record
     * gets emitted when taking the union of variants that contain genotypes. The list must be passed as a
     * comma-separated string listing the names of the variant input files. The list must be complete and include all
     * variant inputs that are being provided to the tool. Use name tags for best results.
     */
    @Argument(fullName="priority_list", shortName="priority", doc="Ordered list specifying priority for merging", optional = true)
    public String PRIORITY_STRING = null;

    @Argument(fullName="printComplexMerges", shortName="printComplexMerges", doc="Emit interesting sites requiring complex compatibility merging to file", optional = true)
    public boolean printComplexMerges = false;

    /**
     * If enabled, this flag causes filtered variants (i.e. variant records where the FILTER field is populated by
     * something other than PASS or a dot) to be omitted from the output.
     */
    @Argument(fullName="filteredAreUncalled", shortName="filteredAreUncalled", doc="Treat filtered variants as uncalled", optional = true)
    public boolean filteredAreUncalled = false;

    /**
     * If this flag is enabled, the INFO, FORMAT and sample-level (genotype) fields will not be emitted to the output file.
     */
    @Argument(fullName="minimalVCF", shortName="minimalVCF", doc="Emit a sites-only file", optional = true)
    public boolean minimalVCF = false;

    /**
     * Exclude sites that do not contain any called ALT alleles in the merged callset. The evaluation is made after the
     * merging procedure is complete.
     */
    @Argument(fullName="excludeNonVariants", shortName="env", doc="Exclude sites where no variation is present after merging", optional = true)
    public boolean EXCLUDE_NON_VARIANTS = false;

    /**
     * Key used in the INFO key=value tag emitted describing which set(s) the combined record came from
     * (e.g. set=control). This provides the option to override the default naming, so instead of set=control you could
     * have it be origin=control, or any other word you want that is not already an INFO field attribute. Set this to
     * 'null' if you don't want the set attribute emitted at all.
     */
    @Argument(fullName="setKey", shortName="setKey", doc="Key name for the set attribute", optional = true)
    public String SET_KEY = "set";

    /**
     * This option allows you to perform a simple merge (concatenation) to combine the VCFs, drastically reducing
     * runtime. Note that in many cases where you think you want to use this option, you may want to check out the
     * CatVariants tool instead, because CatVariants provides the same functionality, but does so even more efficiently.
     */
    @Argument(fullName="assumeIdenticalSamples", shortName="assumeIdenticalSamples", doc="Assume input VCFs have identical sample sets and disjoint calls", optional = true)
    public boolean ASSUME_IDENTICAL_SAMPLES = false;

    /**
     * Sites that are present in fewer than this number of inputs will be ignored. This is a convenient way to build
     * a collection of common variants and exclude rare variants.
     */
    @Argument(fullName="minimumN", shortName="minN", doc="Minimum number of input files the site must be observed in to be included", optional = true)
    public int minimumN = 1;

    /**
     * By default, this tool writes the command line that was used in the header of the output VCF file. This flag
     * enables you to override that behavior . This is most often useful when combining variants for dozens or
     * hundreds of smaller VCFs iteratively, to avoid cluttering the header with a lot of command lines.
     */
    @Argument(fullName="suppressCommandLineHeader", shortName="suppressCommandLineHeader", doc="Do not output the command line to the header", optional = true)
    public boolean SUPPRESS_COMMAND_LINE_HEADER = false;

    /**
     * By default, the INFO field of the merged variant record only contains the INFO field attributes for which all
     * original overlapping records had the same values. Discordant attributes are therefore discarded. This flag allows you to
     * override that behavior and simply copy over the INFO field contents of whichever record had the highest AC value.
     */
    @Argument(fullName="mergeInfoWithMaxAC", shortName="mergeInfoWithMaxAC", doc="Use the INFO content of the record with the highest AC", optional = true)
    public boolean MERGE_INFO_WITH_MAX_AC = false;

    private List<String> priority = null;

    /** Optimization to strip out genotypes before merging if we are doing a sites_only output */
    private boolean sitesOnlyVCF = false;
    private Set<String> samples;

    private VariantContextWriter vcfWriter = null;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        VCFHeader vcfHeader = getHeaderForVariants();

        Utils.nonNull(out);
        IOUtil.assertFileIsWritable(out);

        validateAnnotateUnionArguments();

        Map<String, VCFHeader> headers = new HashMap<>();
        getDrivingVariantsFeatureInputs().forEach(x -> headers.put(x.getName(), (VCFHeader)getHeaderForFeatures(x)));
        final boolean sampleNamesAreUnique = areSamplesUnique(headers.values());

        if (genotypeMergeOption == null && !ASSUME_IDENTICAL_SAMPLES) {
            if (!sampleNamesAreUnique)
                throw new UserException.BadInput("Duplicate sample names were discovered but no genotypemergeoption was supplied. " +
                        "To combine samples without merging, specify --genotypemergeoption UNIQUIFY. Merging duplicate samples " +
                        "without specified priority is unsupported, but can be achieved by specifying --genotypemergeoption UNSORTED.");
            else
                genotypeMergeOption = GATKVariantContextUtils.GenotypeMergeType.UNSORTED;
        }

        if ( PRIORITY_STRING == null && genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE) {
            logger.info("Priority string is not provided, using arbitrary genotyping order: " + priority);
        }

        if (genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE && !sampleNamesAreUnique) {
            throw new IllegalStateException("REQUIRE_UNIQUE sample names is true but duplicate names were discovered.");
        }

        samples = sitesOnlyVCF ? Collections.emptySet() : VcfUtils.getSortedSampleSet(headers, genotypeMergeOption);

        if ( SET_KEY.toLowerCase().equals("null") )
            SET_KEY = null;

        if ( SET_KEY != null ) {
            vcfHeader.addMetaDataLine(new VCFInfoHeaderLine(SET_KEY, 1, VCFHeaderLineType.String, "Source VCF for the merged record in MergeVcfsAndGenotypes"));
        }

        if ( !ASSUME_IDENTICAL_SAMPLES ) {
            vcfHeader.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
            vcfHeader.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
            vcfHeader.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        }

        vcfHeader = new VCFHeader(vcfHeader.getMetaDataInInputOrder(), sitesOnlyVCF ? Collections.emptySet() : samples);

        vcfWriter = createVCFWriter(out);
        vcfHeader.setWriteCommandLine(!SUPPRESS_COMMAND_LINE_HEADER);
        vcfWriter.writeHeader(vcfHeader);
    }

    public static boolean areSamplesUnique(Collection<VCFHeader> headers) {
        Set<String> samples = new HashSet<>();
        for (VCFHeader header : headers) {
            if (CollectionUtils.containsAny(samples, header.getGenotypeSamples())) {
                return false;
            }

            samples.addAll(header.getGenotypeSamples());
        }

        return true;
    }

    private void validateAnnotateUnionArguments() {
        Set<String> rodNames = getDrivingVariantsFeatureInputs().stream().map(FeatureInput::getName).collect(Collectors.toSet());

        if ( genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE && PRIORITY_STRING == null )
            throw new UserException.BadInput("Priority string must be provided if you want to prioritize genotypes");

        if ( PRIORITY_STRING != null){
            priority = new ArrayList<>(Arrays.asList(PRIORITY_STRING.split(",")));
            if ( rodNames.size() != priority.size() )
                throw new UserException.BadInput("The priority list must contain exactly one rod binding per ROD provided to the GATK: rodNames=" + rodNames + " priority=" + priority);

            if ( ! rodNames.containsAll(priority) )
                throw new UserException.BadInput("Not all priority elements provided as input RODs: " + PRIORITY_STRING);
        }
    }

    @Override
    public void apply(List<VariantContext> vcs, ReferenceContext referenceContext, List<ReadsContext> readsContexts) {
        final int totalVCFs = getDrivingVariantsFeatureInputs().size();

        if ( sitesOnlyVCF ) {
            vcs = new ArrayList<>(VariantContextUtils.sitesOnlyVariantContexts(vcs));
        }

        if ( ASSUME_IDENTICAL_SAMPLES ) {
            for ( final VariantContext vc : vcs ) {
                vcfWriter.add(vc);
            }

            return;
        }

        int numFilteredRecords = 0;
        for (final VariantContext vc : vcs) {
            if (vc.filtersWereApplied() && vc.isFiltered())
                numFilteredRecords++;
        }

        if (minimumN > 1 && (vcs.size() - numFilteredRecords < minimumN))
            return;

        final List<VariantContext> mergedVCs = new ArrayList<>();

        if (multipleAllelesMergeType == MultipleAllelesMergeType.BY_TYPE) {
            final Map<VariantContext.Type, List<VariantContext>> VCsByType = separateVariantContextsByType(vcs);

            // TODO -- clean this up in a refactoring
            // merge NO_VARIATION into another type of variant (based on the ordering in VariantContext.Type)
            if ( VCsByType.containsKey(VariantContext.Type.NO_VARIATION) && VCsByType.size() > 1 ) {
                final List<VariantContext> refs = VCsByType.remove(VariantContext.Type.NO_VARIATION);
                for ( final VariantContext.Type type : VariantContext.Type.values() ) {
                    if ( VCsByType.containsKey(type) ) {
                        VCsByType.get(type).addAll(refs);
                        break;
                    }
                }
            }

            // iterate over the types so that it's deterministic
            for (final VariantContext.Type type : VariantContext.Type.values()) {
                // make sure that it is a variant or in case it is not, that we want to include the sites with no variants
                if (!EXCLUDE_NON_VARIANTS || !type.equals(VariantContext.Type.NO_VARIATION)) {
                    if (VCsByType.containsKey(type)) {
                        mergedVCs.add(GATKVariantContextUtils.simpleMerge(VCsByType.get(type), priority, totalVCFs, filteredRecordsMergeType, genotypeMergeOption, filteredAreUncalled));
                    }
                }
            }
        }
        else if (multipleAllelesMergeType == MultipleAllelesMergeType.MIX_TYPES) {
            mergedVCs.add(GATKVariantContextUtils.simpleMerge(vcs, priority, totalVCFs, filteredRecordsMergeType, genotypeMergeOption, filteredAreUncalled));
        }
        else {
            logger.warn("Ignoring all records at site " + referenceContext.getContig() + ":" + referenceContext.getStart());
        }

        for ( final VariantContext mergedVC : mergedVCs ) {
            // only operate at the start of events
            if ( mergedVC == null )
                continue;

            if (hasNonRefSymbolicAllele(mergedVC)) {
                throw new UserException("MergeVcfsAndGenotypes should not be used to merge gVCFs produced by the HaplotypeCaller; use CombineGVCFs instead");
            }

            final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);
            // re-compute chromosome counts
            VariantContextUtils.calculateChromosomeCounts(builder, false);

            if ( minimalVCF )
                pruneVariantContext(builder, Collections.singleton(SET_KEY));
            final VariantContext vc = builder.make();
            if( !EXCLUDE_NON_VARIANTS || vc.isPolymorphicInSamples() )
                vcfWriter.add(builder.make());
        }
    }

    private boolean hasNonRefSymbolicAllele(VariantContext vc) {
        for (Allele allele : vc.getAlleles()) {
            if (allele.isSymbolic() && allele.getBaseString().equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME)) {
                return true;
            }
        }

        return false;
    }

    // Ported from GATK3 GATKVariantContextUtils
    private void pruneVariantContext(final VariantContextBuilder builder, Collection<String> keysToPreserve ) {
        final VariantContext vc = builder.make();
        if ( keysToPreserve == null ) keysToPreserve = Collections.emptyList();

        // VC info
        final Map<String, Object> attributes = subsetAttributes(vc.getCommonInfo(), keysToPreserve);

        // Genotypes
        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype g : vc.getGenotypes() ) {
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            // remove AD, DP, PL, and all extended attributes, keeping just GT and GQ
            gb.noAD().noDP().noPL().noAttributes();
            genotypes.add(gb.make());
        }

        builder.genotypes(genotypes).attributes(attributes);
    }

    // Ported from GATK3 GATKVariantContextUtils
    private Map<String, Object> subsetAttributes(final CommonInfo igc, final Collection<String> keysToPreserve) {
        Map<String, Object> attributes = new HashMap<>(keysToPreserve.size());
        for ( final String key : keysToPreserve  ) {
            if ( igc.hasAttribute(key) )
                attributes.put(key, igc.getAttribute(key));
        }
        return attributes;
    }

    public enum MultipleAllelesMergeType {
        /**
         * Combine only alleles of the same type (SNP, indel, etc.) into a single VCF record.
         */
        BY_TYPE,
        /**
         * Merge all allele types at the same start position into the same VCF record.
         */
        MIX_TYPES
    }

    private Map<VariantContext.Type, List<VariantContext>> separateVariantContextsByType( final Collection<VariantContext> VCs ) {
        if( VCs == null ) { throw new IllegalArgumentException("VCs cannot be null."); }

        final HashMap<VariantContext.Type, List<VariantContext>> mappedVCs = new HashMap<>();
        for ( final VariantContext vc : VCs ) {
            VariantContext.Type vcType = vc.getType();

            // look at previous variant contexts of different type. If:
            // a) otherVC has alleles which are subset of vc, remove otherVC from its list and add otherVC to vc's list
            // b) vc has alleles which are subset of otherVC. Then, add vc to otherVC's type list (rather, do nothing since vc will be added automatically to its list)
            // c) neither: do nothing, just add vc to its own list
            boolean addtoOwnList = true;
            for (final VariantContext.Type type : VariantContext.Type.values()) {
                if (type.equals(vcType))
                    continue;

                if (!mappedVCs.containsKey(type))
                    continue;

                List<VariantContext> vcList = mappedVCs.get(type);
                for (int k=0; k <  vcList.size(); k++) {
                    VariantContext otherVC = vcList.get(k);
                    if (allelesAreSubset(otherVC,vc)) {
                        // otherVC has a type different than vc and its alleles are a subset of vc: remove otherVC from its list and add it to vc's type list
                        vcList.remove(k);
                        // avoid having empty lists
                        if (vcList.isEmpty())
                            mappedVCs.remove(type);
                        if ( !mappedVCs.containsKey(vcType) )
                            mappedVCs.put(vcType, new ArrayList<>());
                        mappedVCs.get(vcType).add(otherVC);
                        break;
                    }
                    else if (allelesAreSubset(vc,otherVC)) {
                        // vc has a type different than otherVC and its alleles are a subset of VC: add vc to otherVC's type list and don't add to its own
                        mappedVCs.get(type).add(vc);
                        addtoOwnList = false;
                        break;
                    }
                }
            }
            if (addtoOwnList) {
                if ( !mappedVCs.containsKey(vcType) )
                    mappedVCs.put(vcType, new ArrayList<VariantContext>());
                mappedVCs.get(vcType).add(vc);
            }
        }

        return mappedVCs;
    }

    private boolean allelesAreSubset(VariantContext vc1, VariantContext vc2) {
        // if all alleles of vc1 are a contained in alleles of vc2, return true
        if (!vc1.getReference().equals(vc2.getReference()))
            return false;

        for (final Allele a :vc1.getAlternateAlleles()) {
            if (!vc2.getAlternateAlleles().contains(a))
                return false;
        }

        return true;
    }

    @Override
    public void closeTool() {
        super.closeTool();

        if (vcfWriter != null)
            vcfWriter.close();
    }
}
