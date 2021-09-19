package com.github.discvrseq.walkers.annotator;

import com.github.discvrseq.tools.VariantManipulationProgramGroup;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKVariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator;

import java.util.*;

/**
 * Annotate variant calls with context information. This is an extension of GATK's VariantAnnotator, and adds a number of custom optional annotations. Usage is essentially identical as GATK's tool.
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to annotate and optionally one or more BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * An annotated VCF.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <br />
 *
 * <h4>Annotate a VCF with dbSNP IDs and depth of coverage for each sample</h4>
 * <pre>
 *   VariantAnnotator \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -A Coverage \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>Annotate a VCF with allele frequency by an external resource. Annotation will only occur if there is allele concordance between the resource and the input VCF </h4>
 * <pre>
 *   VariantAnnotator \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -L anotherInput.vcf \
 *   --resource foo:resource.vcf \
 *   -E foo.AF \
 *   --resource-allele-concordance
 * </pre>
 *
 * Additional annotations include:
 * - GenotypeConcordance: flag genotypes with discordant values relative to a reference VCF
 * - GenotypeConcordanceBySite: Annotate the number of discordant genotypes per site, relative to a reference VCF
 * - MendelianViolationCount: Uses a more extensive check for MVs than default GATK. This annotates the number of MVs and IDs of samples with MVs for each site.
 * - MendelianViolationCountBySample: This is a genotype-level annotation that flags the number of MVs detected per sample.
 *
 */
@CommandLineProgramProperties(summary="Tool for adding annotations to VCF files",
        oneLineSummary = "Tool for adding annotations to VCF files",
        programGroup = VariantManipulationProgramGroup.class)

@DocumentedFeature
public class DiscvrVariantAnnotator extends VariantAnnotator {
    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        return Collections.singletonList(new DiscvrAnnotationPluginDescriptor());
    }

    @Override
    public Collection<Annotation> makeVariantAnnotations() {
        Set<Annotation> ret = new HashSet<>();

        DiscvrAnnotationPluginDescriptor annotationPlugin = getCommandLineParser().getPluginDescriptor(DiscvrAnnotationPluginDescriptor.class);
        ret.addAll(annotationPlugin.getResolvedInstances());

        annotationPlugin.genotypeConcordanceArgumentCollection.featureManager = this.features;

        return ret;
    }

    private static class DiscvrAnnotationPluginDescriptor extends GATKAnnotationPluginDescriptor
    {
        @ArgumentCollection
        public GenotypeConcordanceArgumentCollection genotypeConcordanceArgumentCollection = new GenotypeConcordanceArgumentCollection();

        @ArgumentCollection
        public MendelianViolationArgumentCollection mendelianViolationArgumentCollection = new MendelianViolationArgumentCollection();

        public DiscvrAnnotationPluginDescriptor()
        {
            super(new DefaultGATKVariantAnnotationArgumentCollection(), Collections.emptyList(), Collections.emptyList());
        }

        @Override
        public List<String> getPackageNames() {
            List<String> packageNames = new ArrayList<>(super.getPackageNames());
            packageNames.add(getClass().getPackage().getName());

            return packageNames;
        }

        @Override
        public void validateAndResolvePlugins() throws CommandLineException {
            super.validateAndResolvePlugins();

            getResolvedInstances().stream()
                    .filter(GenotypeConcordanceArgumentCollection.UsesGenotypeConcordanceArgumentCollection.class::isInstance)
                    .map(a -> (GenotypeConcordanceArgumentCollection.UsesGenotypeConcordanceArgumentCollection) a)
                    .forEach(a -> {
                        a.setArgumentCollection(genotypeConcordanceArgumentCollection);
                        genotypeConcordanceArgumentCollection.validateArguments();
                    });

            getResolvedInstances().stream()
                    .filter(MendelianViolationArgumentCollection.UsesMendelianViolationArgumentCollection.class::isInstance)
                    .map(a -> (MendelianViolationArgumentCollection.UsesMendelianViolationArgumentCollection) a)
                    .forEach(a -> {
                        a.setArgumentCollection(mendelianViolationArgumentCollection);
                        mendelianViolationArgumentCollection.validateArguments();
                    });
        }
    }
}
