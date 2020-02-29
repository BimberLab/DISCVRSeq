package com.github.discvrseq.util;

import org.apache.commons.text.similarity.LevenshteinDistance;

public class SequenceMatcher {
    private static final LevenshteinDistance levenshteinDistance = LevenshteinDistance.getDefaultInstance();

    public static Integer fuzzyMatch(String query, String readSeq, int editDistance) {
        int windows = readSeq.length() - query.length();

        int i = 0;
        while (i < windows) {
            CharSequence test = readSeq.subSequence(i, i + query.length());
            if (levenshteinDistance.apply(query, test) <= editDistance) {
                return i;
            }

            i++;
        }

        return null;
    }
}
