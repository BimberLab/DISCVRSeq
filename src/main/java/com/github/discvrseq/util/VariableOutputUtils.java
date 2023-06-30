package com.github.discvrseq.util;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.lang.reflect.Array;
import java.util.*;
import java.util.function.Function;

public class VariableOutputUtils {
    private static final String MISSING_DATA = "NA";

    public static List<List<String>> extractFields(final VariantContext vc, final VCFHeader inputHeader, final List<String> fieldsToTake, final List<String> asFieldsToTake, final boolean errorIfMissingData, final boolean splitMultiAllelic) {
        final int numRecordsToProduce = splitMultiAllelic ? vc.getAlternateAlleles().size() : 1;
        final List<List<String>> records = new ArrayList<>(numRecordsToProduce);

        for ( int i = 0; i < numRecordsToProduce; i++ ) {
            records.add(new ArrayList<>());
        }

        for ( final String field : fieldsToTake ) {
            if ( splitMultiAllelic && field.equals("ALT") ) { // we need to special case the ALT field when splitting out multi-allelic records
                addFieldValue(splitAltAlleles(vc), records);
            } else if ( getters.containsKey(field) ) {
                addFieldValue(getters.get(field).apply(vc), records);
            } else if ( vc.hasAttribute(field) ) {
                addFieldValue(vc.getAttribute(field, null), records);
            } else if ( isWildCard(field) ) {
                final SortedSet<String> wildVals = new TreeSet<>();
                for ( final Map.Entry<String,Object> elt : vc.getAttributes().entrySet()) {
                    if ( elt.getKey().startsWith(field.substring(0, field.length() - 1)) ) {
                        wildVals.add(elt.getValue().toString());
                    }
                }

                final String val = wildVals.isEmpty() ? MISSING_DATA : Utils.join(",", wildVals);

                addFieldValue(val, records);
            } else {
                handleMissingData(errorIfMissingData, field, records, vc);
            }
        }

        for ( final String field : asFieldsToTake) {
            if (vc.hasAttribute(field)) {
                if (splitMultiAllelic) {
                    addAlleleSpecificFieldValue(Arrays.asList(vc.getAttributeAsString(field, ".").replace("[", "").replace("]", "").split(",")), records, inputHeader.getInfoHeaderLine(field).getCountType());

                } else {
                    addFieldValue(vc.getAttributeAsString(field, ".").replace("[","").replace("]",""), records);
                }
            } else {
                handleMissingData(errorIfMissingData, field, records, vc);
            }
        }

        return records;
    }

    private static void addFieldValue(final Object val, final List<List<String>> result) {
        final int numResultRecords = result.size();

        // if we're trying to create a single output record, add it
        if ( numResultRecords == 1 ) {
            result.get(0).add(prettyPrintObject(val));
        }
        // if this field is a list of the proper size, add the appropriate entry to each record
        else if ( (val instanceof List) && ((List)val).size() == numResultRecords ) {
            final List<?> list = (List<?>)val;
            for ( int i = 0; i < numResultRecords; i++ ) {
                result.get(i).add(list.get(i).toString());
            }
        }
        // otherwise, add the original value to all of the records
        else {
            final String valStr = prettyPrintObject(val);
            for ( final List<String> record : result ) {
                record.add(valStr);
            }
        }
    }

    private static void addAlleleSpecificFieldValue(final Object val, final List<List<String>> result, final VCFHeaderLineCount alleleCount) {
        if (val instanceof List && alleleCount.equals(VCFHeaderLineCount.R)) {
            final List<?> myList = (List<?>) val;
            addFieldValue(new ArrayList<>(myList.subList(1, myList.size())), result);
        }
        else {
            addFieldValue(val, result);
        }
    }

    private static String prettyPrintObject(final Object val) {
        if ( val == null ) {
            return "";
        }

        if ( val instanceof List ) {
            return prettyPrintObject(((List) val).toArray());
        }

        if ( !val.getClass().isArray() ) {
            return val.toString();
        }

        final int length = Array.getLength(val);
        if ( length == 0 ) {
            return "";
        }

        final StringBuilder sb = new StringBuilder(prettyPrintObject(Array.get(val, 0)));
        for ( int i = 1; i < length; i++ ) {
            sb.append(",");
            sb.append(prettyPrintObject(Array.get(val, i)));
        }
        return sb.toString();
    }

    public static boolean hasGetter(String key) {
        return getters.containsKey(key);
    }

    private static final Map<String, Function<VariantContext, String>> getters = new LinkedHashMap<>();

    static {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", vc -> vc.getContig());
        getters.put("POS", vc -> Integer.toString(vc.getStart()));
        getters.put("REF", vc -> vc.getReference().getDisplayString());
        getters.put("ALT", vc -> {
            final StringBuilder x = new StringBuilder();
            final int n = vc.getAlternateAlleles().size();
            if ( n == 0 ) {
                return ".";
            }

            for ( int i = 0; i < n; i++ ) {
                if ( i != 0 ) {
                    x.append(",");
                }
                x.append(vc.getAlternateAllele(i));
            }
            return x.toString();
        });
        getters.put("EVENTLENGTH", vc -> {
            int maxLength = 0;
            for ( final Allele a : vc.getAlternateAlleles() ) {
                final int length = a.length() - vc.getReference().length();
                if( Math.abs(length) > Math.abs(maxLength) ) { maxLength = length; }
            }
            return Integer.toString(maxLength);
        });
        getters.put("QUAL", vc -> Double.toString(vc.getPhredScaledQual()));
        getters.put("TRANSITION", vc -> {
            if ( vc.isSNP() && vc.isBiallelic() ) {
                return GATKVariantContextUtils.isTransition(vc) ? "1" : "0";
            } else {
                return "-1";
            }
        });
        getters.put("FILTER", vc -> vc.isNotFiltered() ? "PASS" : Utils.join(",", vc.getFilters()));
        getters.put("ID", vc -> vc.getID());
        getters.put("HET", vc -> Integer.toString(vc.getHetCount()));
        getters.put("HOM-REF", vc -> Integer.toString(vc.getHomRefCount()));
        getters.put("HOM-VAR", vc -> Integer.toString(vc.getHomVarCount()));
        getters.put("NO-CALL", vc -> Integer.toString(vc.getNoCallCount()));
        getters.put("TYPE", vc -> vc.getType().toString());
        getters.put("VAR", vc -> Integer.toString(vc.getHetCount() + vc.getHomVarCount()));
        getters.put("NSAMPLES", vc -> Integer.toString(vc.getNSamples()));
        getters.put("NCALLED", vc -> Integer.toString(vc.getNSamples() - vc.getNoCallCount()));
        getters.put("MULTI-ALLELIC", vc -> Boolean.toString(vc.getAlternateAlleles().size() > 1));
        getters.put("SAMPLE_NAME", vc -> vc.getGenotype(0).getSampleName());
    }

    private static Object splitAltAlleles(final VariantContext vc) {
        final int numAltAlleles = vc.getAlternateAlleles().size();
        if ( numAltAlleles == 1 ) {
            return vc.getAlternateAllele(0);
        }

        return vc.getAlternateAlleles();
    }

    private static boolean isWildCard(final String s) {
        return s.endsWith("*");
    }

    private static void handleMissingData(final boolean errorIfMissingData, final String field, final List<List<String>> records, final VariantContext vc) {
        if (errorIfMissingData) {
            throw new UserException(String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc));
        } else {
            addFieldValue(MISSING_DATA, records);
        }
    }
}
