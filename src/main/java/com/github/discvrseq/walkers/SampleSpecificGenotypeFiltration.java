package com.github.discvrseq.walkers;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.util.*;

/**
 * This tool is a specialized version of VariantFiltration. It takes a tab-delimited text file mapping sample name to group, and then allows the user to perform group-specific JEXL genotype filtration.
 * An example use-case would be a VCF with a mixture of WGS and WXS data, where samples from each technology require different filter thresholds.
 *
 * <h3>Usage example:</h3>
 * <pre>
 *  java -jar DISCVRseq.jar SampleSpecificGenotypeFiltration \
 *     -R currentGenome.fasta \
 *     -V myVCF.vcf \
 *     --sample-map sampleMap.txt \
 *     -O output.vcf.gz \
 *     --genotype-filter-expression 'Set1:DP < 10'
 *     --genotype-filter-name 'DP-LT10'
 *     --genotype-filter-expression 'Set2:DP < 10'
 *     --genotype-filter-name 'DP-LT10'
 *     --genotype-filter-expression 'Set2:DP > 50'
 *     --genotype-filter-name 'DP-GT50'
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "This tool is a specialized version of VariantFiltration. It takes a text file mapping sample name to group, and then allows the user to perform group-specific JEXL genotype filtration.",
        oneLineSummary = "Allows filtration of genotypes with sample-specific expressions",
        programGroup = VariantManipulationProgramGroup.class
)
public class SampleSpecificGenotypeFiltration extends VariantWalker {
    private static final String FILTER_DELIMITER = ";";

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public GATKPath out = null;

    @Argument(doc="File with sample to setName mapping", fullName = "sample-map", shortName = "sm", optional = false)
    public GATKPath sampleMap = null;

    @Argument(fullName= VariantFiltration.GENOTYPE_FILTER_EXPRESSION_LONG_NAME, shortName="G-filter", doc="One or more expressions used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info). These must be in the format \"SetName:<JEXL>\", where SetName matches a sample group defined by --sampleMap", optional=false)
    public List<String> genotypeFilterExpressions = new ArrayList<>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     */
    @Argument(fullName=VariantFiltration.GENOTYPE_FILTER_NAME_LONG_NAME, shortName="G-filter-name", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", optional=false)
    public List<String> genotypeFilterNames = new ArrayList<>();

    @Argument(fullName=VariantFiltration.MISSING_VAL_LONG_NAME, doc="When evaluating the JEXL expressions, missing values should be considered failing the expression", optional=true)
    public Boolean failMissingValues = false;

    private Map<String, String> sampleToSetName;
    private final Map<String, List<VariantContextUtils.JexlVCMatchExp>> setToFilter = new HashMap<>();
    private VariantContextWriter writer;

    private JexlMissingValueTreatment howToTreatMissingValues;

    @Override
    public void onTraversalStart() {
        sampleToSetName = loadSampleNameMapFile();

        if (genotypeFilterExpressions.size() != genotypeFilterNames.size()) {
            throw new GATKException("The same number of genotype filters and filter names must be provided");
        }

        // Initialize per-set filters
        Map<String, Pair<List<String>, List<String>>> setToFilterStrings = new HashMap<>();
        Set<String> allSetNames = new HashSet<>(sampleToSetName.values());
        for (int i=0;i<genotypeFilterExpressions.size();i++) {
            String expr = genotypeFilterExpressions.get(i);
            String[] tokens = expr.split(":");

            String setName = tokens[0];
            String genotypeExpr = StringUtils.join(Arrays.copyOfRange(tokens, 1, tokens.length));

            if (!allSetNames.contains(setName)) {
                throw new GATKException("Unknown set name for filter: " + expr);
            }

            Pair<List<String>, List<String>> filters = setToFilterStrings.get(setName);
            if (filters == null) {
                filters = Pair.of(new ArrayList<>(), new ArrayList<>());
            }

            filters.getLeft().add(genotypeFilterNames.get(i));
            filters.getRight().add(genotypeExpr);

            setToFilterStrings.put(setName, filters);
        }

        setToFilterStrings.forEach((setName, filterList) -> {
            setToFilter.put(setName, VariantContextUtils.initializeMatchExps(filterList.getLeft(), filterList.getRight()));
        });

        howToTreatMissingValues = failMissingValues ? JexlMissingValueTreatment.TREAT_AS_MATCH : JexlMissingValueTreatment.TREAT_AS_MISMATCH;

        initializeVcfWriter();
    }

    private void initializeVcfWriter() {
        writer = createVCFWriter(out);

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new LinkedHashSet<>(getHeaderForVariants().getMetaDataInInputOrder());

        // need AC, AN and AF since output if set filtered genotypes to no-call
        // If setting filtered genotypes to no-call, then allele counts (AC, AN and AF ) will be recomputed and these annotations
        // need to be included in the header
        GATKVariantContextUtils.addChromosomeCountsToHeader(hInfo);

        hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));

        try {
            setToFilter.values().stream().flatMap(List::stream).forEach(exp -> {
                hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
            });
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput(e.getMessage());
        }

        writer.writeHeader(new VCFHeader(hInfo, getHeaderForVariants().getGenotypeSamples()));
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder builder = new VariantContextBuilder(variant);

        // make new Genotypes based on filters
        for (String setName : setToFilter.keySet()) {
            List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps = setToFilter.get(setName);
            if ( !genotypeFilterExps.isEmpty()) {
                GATKVariantContextUtils.setFilteredGenotypeToNocall(builder, variant, true, this::getGenotypeFilters);
            }
        }

        writer.add(builder.make());
    }

    private List<String> getGenotypeFilters(final VariantContext vc, final Genotype g) {
        final List<String> filters = new ArrayList<>();
        if (g.isFiltered()) {
            filters.addAll(Arrays.asList(g.getFilters().split(FILTER_DELIMITER)));
        }

        // Add if expression filters the variant context
        String setName = sampleToSetName.get(g.getSampleName());
        if (setName == null) {
            return filters;
        }

        List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps = setToFilter.get(setName);
        if (genotypeFilterExps == null || genotypeFilterExps.isEmpty()) {
            return filters;
        }

        for (final VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps) {
            if (matchesFilter(vc, g, exp)) {
                filters.add(exp.name);
            }
        }

        return filters;
    }

    private boolean matchesFilter(final VariantContext vc, final Genotype g, final VariantContextUtils.JexlVCMatchExp exp) {
        return VariantContextUtils.match(vc, g, exp, howToTreatMissingValues);
    }

    private LinkedHashMap<String, String> loadSampleNameMapFile() {
        IOUtil.assertFileIsReadable(sampleMap.toPath());

        try {
            final List<String> lines = Files.readAllLines(sampleMap.toPath());
            if (lines.isEmpty()) {
                throw new UserException.BadInput( "At least 1 sample is required but none were found in the sample mapping file");
            }

            final LinkedHashMap<String, String> sampleToSetName = new LinkedHashMap<>();
            for ( final String line : lines) {
                final String[] split = line.split("\\t",-1);
                if (split.length != 2) {
                    throw new UserException.BadInput("Expected a file with 2 fields per line in the format\nSample\tSetName\n but found line: \""
                            + line +"\" with "+split.length+" fields");
                }
                if ( !split[0].trim().equals(split[0]) || split[0].trim().isEmpty()
                        || split[1].trim().isEmpty()) {
                    throw new UserException.BadInput("Expected a file of format\nSample\tSetName\n but found line: '" + line + "'\nValid sample names must be non-empty strings that cannot begin or end with whitespace and valid set names must be non-empty and not all whitespace");
                }

                final String sample = split[0];
                final String setName = split[1].trim();

                sampleToSetName.put(sample, setName);
            }
            return sampleToSetName;
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(sampleMap, "exception while reading sample->setName mapping file",  e);
        }
    }

    @Override
    public void closeTool() {
        if (writer != null){
            writer.close();
        }
    }
}
