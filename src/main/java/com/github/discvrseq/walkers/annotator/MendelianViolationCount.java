package com.github.discvrseq.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.PedigreeAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.Sex;

import java.util.*;

/**
 * Annotate the number of subjects with a Mendelian Violation
 *
 * <p>This annotation uses the likelihoods of the genotype calls to assess whether a site is transmitted from parents to offspring according to Mendelian rules. The output is the likelihood of the site being a Mendelian violation, which can be tentatively interpreted either as an indication of error (in the genotype calls) or as a possible <em><de novo/em> mutation. The higher the output value, the more likely there is to be a Mendelian violation. Note that only positive values indicating likely MVs will be annotated; if the value for a given site is negative (indicating that there is no violation) the annotation is not written to the file.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>The defaults to a hard cutoff of 10 for the phred scale genotype quality.  If a given sample is below this threshold, the genotype is assumed to be no call.</li>
 *     <li>This annotation requires a valid pedigree file.</li>
 * </ul>
 *
 */
public class MendelianViolationCount extends PedigreeAnnotation implements InfoFieldAnnotation {
    public static final String MV_NUM = "MV_NUM";
    public static final String MV_SAMPLES = "MV_SAMPLES";
    private SampleDB sampleDB = null;

    @ArgumentCollection
    public MendelianViolationArgumentCollection args = new MendelianViolationArgumentCollection();

    public MendelianViolationCount()
    {
        super((Set<String>) null);
    }

    public MendelianViolationCount(final GATKPath pedigreeFile){
        super(pedigreeFile);
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext referenceContext, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods) {
        Map<String,Object> attributeMap = new HashMap<>();
        if (sampleDB != null) {
            int totalViolations = 0;
            Set<String> violations = new HashSet<>();
            for (Sample s : sampleDB.getSamples()) {
                int count = countViolations(s, vc, args.minGenotypeQuality);
                totalViolations += count;
                if (count > 0) {
                    violations.add(s.getID());
                }
            }

            attributeMap.put(MV_NUM, totalViolations);
            attributeMap.put(MV_SAMPLES, StringUtils.join(violations.toArray(new String[0]), ","));
        }

        return attributeMap;
    }

    public static int countViolations(Sample subject, VariantContext vc, double minGenotypeQuality) {
        MV ret = getMendelianViolation(subject, vc, minGenotypeQuality);
        if (ret == null){
            return 0;
        }

        return ret.isViolation() ? 1 : 0;
    }

    public static MV getMendelianViolation(Sample subject, VariantContext vc, double minGenotypeQuality) {
        Genotype gChild = vc.getGenotype(subject.getID());
        if (gChild == null || !gChild.isCalled()){
            return null;  //cant make call
        }

        //Count lowQual. Note that if min quality is set to 0, even values with no quality associated are returned
        if (minGenotypeQuality > -1 && gChild.getGQ() < minGenotypeQuality) {
            return null; //cannot make determination
        }

        //until we have improved calling of sex chromosomes, skip this situation
        //chrY is overwhelmingly no-call for females, so we can proceed w/ the check here
        boolean isChrX = vc.getContig().equalsIgnoreCase("X") || vc.getContig().equalsIgnoreCase("chrX");
        if (isChrX && subject.getSex().equals(Sex.MALE)){
            return null;
        }

        Genotype gMom = vc.getGenotype(subject.getMaternalID());
        if (gMom == null)
        {
            gMom = new NoCallGenotype(subject.getMaternalID());
        }
        Genotype gDad = vc.getGenotype(subject.getPaternalID());
        if (gDad == null)
        {
            gDad = new NoCallGenotype(subject.getPaternalID());
        }

        //If the family is all homref, not too interesting
        if (gMom.isHomRef() && gDad.isHomRef() && gChild.isHomRef()) {
            return null;
        }
        else if (!gMom.isCalled() && !gDad.isCalled()) {
            return null;
        }

        return getViolation(gMom, gDad, gChild, minGenotypeQuality);
    }

    private static MV getViolation(final Genotype gMom, final Genotype gDad, final Genotype gChild, double minGenotypeQuality) {
        MV ret = new MV();

        if (gDad.isCalled()){
            ret.fatherIsViolation = testParent(gChild, gDad, minGenotypeQuality);
        }

        if (gMom.isCalled()){
            ret.motherIsViolation = testParent(gChild, gMom, minGenotypeQuality);
        }

        //Both parents have genotype information
        if (gMom.isCalled() && gDad.isCalled()){
            ret.violationCombined = testParents(gChild, gDad, gMom);
        }

        return ret;
    }

    public static class MV {
        public boolean motherIsViolation = false;
        public boolean fatherIsViolation = false;
        public boolean violationCombined = false;

        public boolean isViolation(){
            return motherIsViolation || fatherIsViolation || violationCombined;
        }
    }

    private static boolean testParents(Genotype gChild, Genotype gDad, Genotype gMom){
        return !(gMom.getAlleles().contains(gChild.getAlleles().get(0)) && gDad.getAlleles().contains(gChild.getAlleles().get(1)) || gMom.getAlleles().contains(gChild.getAlleles().get(1)) && gDad.getAlleles().contains(gChild.getAlleles().get(0)));
    }

    private static boolean testParent(Genotype gChild, Genotype gParent, double minGenotypeQuality){
        if (gParent.getGQ() < minGenotypeQuality) {
            return false;
        }

        return (gParent.isHomRef() && gChild.isHomVar()) || (gParent.isHomVar() && gChild.isHomRef()) || (!gParent.getAlleles().contains(gChild.getAllele(0)) && !gParent.getAlleles().contains(gChild.getAllele(1)));
    }

    // return the descriptions used for the VCF INFO meta field
    public List<String> getKeyNames() { return Arrays.asList(MV_NUM, MV_SAMPLES); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(
            new VCFInfoHeaderLine(MV_NUM, 1, VCFHeaderLineType.Integer, "Number of mendelian violations across all samples."),
            new VCFInfoHeaderLine(MV_SAMPLES, 1, VCFHeaderLineType.String, "Samples where a mendelian violation was observed.")
    ); }

    public static class NoCallGenotype extends Genotype {
        private static final long serialVersionUID = 1L;

        private Genotype _orig = null;
        private List<Allele> _alleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

        public NoCallGenotype(String sampleName)
        {
            super(sampleName, ".");
        }

        @Override
        public List<Allele> getAlleles()
        {
            return _alleles;
        }

        @Override
        public Allele getAllele(int i)
        {
            return _alleles.get(i);
        }

        @Override
        public boolean isPhased()
        {
            return false;
        }

        @Override
        public int getDP()
        {
            return _orig == null ? 0 : _orig.getDP();
        }

        @Override
        public int[] getAD()
        {
            return _orig == null ? new int[0] : _orig.getAD();
        }

        @Override
        public int getGQ()
        {
            return _orig == null ? 50 : _orig.getGQ();
        }

        @Override
        public int[] getPL()
        {
            return _orig == null ? new int[0] : _orig.getPL();
        }

        @Override
        public Map<String, Object> getExtendedAttributes()
        {
            return _orig == null ? null : _orig.getExtendedAttributes();
        }
    }
}
