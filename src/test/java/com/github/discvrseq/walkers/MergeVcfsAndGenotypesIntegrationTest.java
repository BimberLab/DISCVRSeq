package com.github.discvrseq.walkers;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class MergeVcfsAndGenotypesIntegrationTest extends BaseIntegrationTest {

    @Test
    public void basicTest() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("--variant:a");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-V:b");
        File input2 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-priority");
        args.addRaw("a,b");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        String fn = "basicTestOutput.vcf";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest("basicTest", this);
    }

    @Test
    public void basicTest2() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("--variant:a");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-V:b");
        File input2 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-V:c");
        File input3 = new File(testBaseDir, "mergeVcf2.vcf");
        ensureVcfIndex(input3);
        args.addRaw(normalizePath(input3));

        args.add("genotypeMergeOption", "PRIORITIZE");

        args.addRaw("-priority");
        args.addRaw("a,b,c");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        String fn = "basicTest2Output.vcf";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        spec.executeTest("basicTest", this);
    }

    @Test
    public void basicTest3() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("R", normalizePath(getHg19Micro()));

        args.addRaw("--variant:a");
        File input = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input);

        args.addRaw(normalizePath(input));

        args.addRaw("-V:b");
        File input2 = new File(testBaseDir, "mergeVcf1.vcf");
        ensureVcfIndex(input2);
        args.addRaw(normalizePath(input2));

        args.addRaw("-priority");
        args.addRaw("a,b");

        args.addRaw("-O");
        args.addRaw("%s");
        args.addRaw("--tmp-dir");
        args.addRaw(getTmpDir());

        IntegrationTestSpec spec = new IntegrationTestSpec(
                args.getString(),
                1, UserException.BadInput.class);

        spec.executeTest("basicTest", this);


        args.add("genotypeMergeOption", "PRIORITIZE");

        String fn = "basicTest3Output.vcf";
        IntegrationTestSpec spec2 = new IntegrationTestSpec(
                args.getString(),
                Arrays.asList(getTestFile(fn).getPath()));

        executeTest(spec, "basicTest");
    }

    private void executeTest(IntegrationTestSpec spec, String name) throws IOException {
        try {
            spec.executeTest(name, this);
        }
        catch (AssertionError e) {
            throw e;
        }
    }

//    public void test1InOut(String file, String md5, String args) {
//         WalkerTestSpec spec = new WalkerTestSpec(
//                 baseTestString(" -priority v1 -V:v1 " + validationDataLocation + file + args),
//                 1,
//                 Arrays.asList(md5));
//         cvExecuteTest("testInOut1--" + file, spec, true);
//    }
//
//    public void combine2(String file1, String file2, String args, String md5) {
//        combine2(file1, file2, args, md5, true);
//    }
//
//    public void combine2(String file1, String file2, String args, String md5, final boolean parallel) {
//         WalkerTestSpec spec = new WalkerTestSpec(
//                 baseTestString(" -priority v1,v2 -V:v1 " + validationDataLocation + file1 + " -V:v2 "+ validationDataLocation + file2 + args),
//                 1,
//                 Arrays.asList(md5));
//         cvExecuteTest("combine2 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec, parallel);
//    }
//
//    public void combineSites(String args, String md5) {
//        String file1 = "1000G_omni2.5.b37.sites.vcf";
//        String file2 = "hapmap_3.3.b37.sites.vcf";
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineVariants --no_cmdline_in_header -o %s -R " + b37KGReference
//                        + " -L 1:1-10,000,000 -V:omni " + validationDataLocation + file1
//                        + " -V:hm3 " + validationDataLocation + file2 + args,
//                1,
//                Arrays.asList(md5));
//        cvExecuteTest("combineSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec, true);
//    }
//
//    public void combinePLs(String file1, String file2, String md5) {
//         WalkerTestSpec spec = new WalkerTestSpec(
//                 "-T CombineVariants --no_cmdline_in_header -o %s -R " + b36KGReference + " -priority v1,v2 -V:v1 " + privateTestDir + file1 + " -V:v2 " + privateTestDir + file2,
//                 1,
//                 Arrays.asList(md5));
//         cvExecuteTest("combine PLs 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec, true);
//    }
//
//    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "e3dbdfa14aefb2f6bd1213287d34a2e5", " -U LENIENT_VCF_PROCESSING"); }
//    @Test public void test2SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "d727fab83b4265859c4a902f6e66ac3d", " -setKey foo -U LENIENT_VCF_PROCESSING"); }
//    @Test public void test3SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "42fc3d2c68415a61ff15e594a63d9349", " -setKey null -U LENIENT_VCF_PROCESSING"); }
//    @Test public void testOfficialCEUPilotCalls() { test1InOut("CEU.trio.2010_03.genotypes.vcf.gz", "a3994d6145bb3813950939238db4c592"); } // official project VCF files in tabix format
//
//    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "e7fd959312e2aff0b4231963ee690aec"); }
//    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "23439a1f0108b57a14e18efe9482cc88"); }
//
//    @Test public void combineWithPLs() { combinePLs("combine.3.vcf", "combine.4.vcf", "27aa46cdb022be3959e7240a0d7ac794"); }
//
//    @Test public void combineTrioCalls() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", "", "9bdda937754e1407183406808f560723"); } // official project VCF files in tabix format
//    @Test public void combineTrioCallsMin() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", " -minimalVCF", "6344953a82a422115bd647ec1d696b94"); } // official project VCF files in tabix format
//    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "51cf4543e46c8434e32c426f59507d1e"); }
//
//    @Test public void combineSNPsAndIndels() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "f9d1d7e6246f0ce9e493357d5b320323"); }
//
//    @Test public void uniqueSNPs() {
//        // parallelism must be disabled because the input VCF is malformed (DB=0) and parallelism actually fixes this which breaks the md5s
//        //both of these files have the YRI trio and merging of duplicate samples without priority must be specified with UNSORTED merge type
//        combine2("pilot2.snps.vcf4.genotypes.vcf", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", " -genotypeMergeOption UNSORTED", "5aece78046bfb7d6ee8dc4d551542e3a", true);
//    }
//
//    @Test public void omniHM3Union() { combineSites(" -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED", "0897efcc0046bd94760315838d4d0fa5"); }
//    @Test public void omniHM3Intersect() { combineSites(" -filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED", "8b12b09a6ec4e3fde2352bbf82637f1e"); }
//
//    @Test public void threeWayWithRefs() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:NA19240_BGI "+validationDataLocation+"NA19240.BGI.RG.vcf" +
//                        " -V:NA19240_ILLUMINA "+validationDataLocation+"NA19240.ILLUMINA.RG.vcf" +
//                        " -V:NA19240_WUGSC "+validationDataLocation+"NA19240.WUGSC.RG.vcf" +
//                        " -V:denovoInfo "+validationDataLocation+"yri_merged_validation_data_240610.annotated.b36.vcf" +
//                        " -setKey centerSet" +
//                        " -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED" +
//                        " -U LENIENT_VCF_PROCESSING" +
//                        " -priority NA19240_BGI,NA19240_ILLUMINA,NA19240_WUGSC,denovoInfo" +
//                        " -genotypeMergeOption UNIQUIFY -L 1"),
//                1,
//                Arrays.asList("8b75e835ed19c06c358a2185cd0e14db"));
//        cvExecuteTest("threeWayWithRefs", spec, true);
//    }
//
//    // complex examples with filtering, indels, and multiple alleles
//    public void combineComplexSites(String args, String md5) {
//        String file1 = "combine.1.vcf";
//        String file2 = "combine.2.vcf";
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineVariants --no_cmdline_in_header -o %s -R " + b37KGReference
//                        + " -V:one " + privateTestDir + file1
//                        + " -V:two " + privateTestDir + file2 + args,
//                1,
//                Arrays.asList(md5));
//        cvExecuteTest("combineComplexSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec, true);
//    }
//
//    @Test public void complexTestFull() { combineComplexSites("", "80c3b7ba39c8a3f3511fc1ea61ecd4da"); }
//    @Test public void complexTestMinimal() { combineComplexSites(" -minimalVCF", "e7da95bbcf3890a4debbfa07cbd646e5"); }
//    @Test public void complexTestSitesOnly() { combineComplexSites(" -sites_only", "d7ee8da6ceee4dd1212122c2e9cab2a6"); }
//    @Test public void complexTestSitesOnlyMinimal() { combineComplexSites(" -sites_only -minimalVCF", "d7ee8da6ceee4dd1212122c2e9cab2a6"); }
//
//    @Test
//    public void combineDBSNPDuplicateSites() {
//         WalkerTestSpec spec = new WalkerTestSpec(
//                 "-T CombineVariants --no_cmdline_in_header -L 1:902000-903000 -o %s -R " + b37KGReference + " -V:v1 " + b37dbSNP132,
//                 1,
//                 Arrays.asList("b0d4b86702b44fc4faa527c34adf6239"));
//         cvExecuteTest("combineDBSNPDuplicateSites:", spec, true);
//    }
//
//    @Test
//    public void combineLeavesUnfilteredRecordsUnfiltered() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineVariants --no_cmdline_in_header -o %s "
//                        + " -R " + b37KGReference
//                        + " -V " + privateTestDir + "combineVariantsLeavesRecordsUnfiltered.vcf",
//                1,
//                Arrays.asList("0f221847e76521250de1abcba535e49c"));
//        cvExecuteTest("combineLeavesUnfilteredRecordsUnfiltered: ", spec, false);
//    }
//
//    @Test
//    public void combiningGVCFsFails() {
//        try {
//            WalkerTestSpec spec = new WalkerTestSpec(
//                    "-T CombineVariants --no_cmdline_in_header -o %s "
//                            + " -R " + b37KGReference
//                            + " -V " + privateTestDir + "gvcfExample1.vcf",
//                    1,
//                    Arrays.asList("FAILFAILFAILFAILFAILFAILFAILFAIL"));
//            executeTest("combiningGVCFsFails", spec);
//        } catch (Exception e) { } // do nothing
//    }
//
//    @Test
//    public void combineSymbolicVariants() {
//        // Just checking that this does not fail, hence no output files and MD5
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineVariants --no_cmdline_in_header -o %s "
//                        + " -R " + hg19ReferenceWithChrPrefixInChromosomeNames
//                        + " -V " + privateTestDir + "WES-chr1.DEL.vcf"
//                        + " -V " + privateTestDir + "WGS-chr1.DEL.vcf"
//                        + " -genotypeMergeOption UNIQUIFY",
//                0,
//                Arrays.asList(""));
//        executeTest("combineSymbolicVariants: ", spec);
//    }
//
//    @Test
//    public void combineSpanningDels() {
//        // Just checking that this does not fail, hence no output files and MD5
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T CombineVariants --no_cmdline_in_header -o %s "
//                        + " -R " + b37KGReference
//                        + " -V " + privateTestDir + "test.spanningdel.combine.1.vcf "
//                        + " -V " + privateTestDir + "test.spanningdel.combine.2.vcf "
//                        + " -genotypeMergeOption UNIQUIFY",
//                0,
//                Arrays.asList(""));
//        executeTest("combineSpanningDels: ", spec);
//    }
}
