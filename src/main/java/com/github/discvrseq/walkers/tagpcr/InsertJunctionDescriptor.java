package com.github.discvrseq.walkers.tagpcr;

import com.github.discvrseq.util.SequenceMatcher;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class InsertJunctionDescriptor {
    private String name;
    private List<String> searchStrings;
    private List<String> searchStringsRC = null;

    private int mismatchesAllowed = 2;
    private boolean invertHitOrientation = false;

    public InsertJunctionDescriptor() {

    }

    public void setMismatchesAllowed(int mismatchesAllowed) {
        this.mismatchesAllowed = mismatchesAllowed;
    }

    private List<String> getSearchStrings(boolean reverseComplement) {
        if (reverseComplement) {
            if (searchStringsRC == null) {
                searchStringsRC = new ArrayList<>();
                searchStrings.forEach(x -> {
                    searchStringsRC.add(SequenceUtil.reverseComplement(x));
                });
            }

            return searchStringsRC;
        }

        return searchStrings;
    }

    public Map<String, IntegrationSiteMapper.JunctionMatch> getMatches(SAMRecord rec, InsertDescriptor id, Logger log, int maxRecordsToStore) {
        Map<String, IntegrationSiteMapper.JunctionMatch> matches = new HashMap<>();

        //Forward orientation.
        for (String query : getSearchStrings(false)) {
            Integer match0 = SequenceMatcher.fuzzyMatch(query, rec.getReadString().toUpperCase(), mismatchesAllowed);
            if (match0 != null) {
                //this is equal to one base after the 1-based end
                final int start = match0 + query.length() + 2;
                int pos = rec.getReferencePositionAtReadPosition(start);

                int i = 1;
                while (pos == 0 && (start + i) <= rec.getReadLength()) {
                    pos = rec.getReferencePositionAtReadPosition(start + i);
                    i++;
                }

                if (pos > 0) {
                    IntegrationSiteMapper.JunctionMatch jm = new IntegrationSiteMapper.JunctionMatch(log, id, this, rec.getContig(), pos, false, rec, maxRecordsToStore);
                    matches.put(jm.getKey(), jm);
                }
            }
        }

        //Reverse:
        for (String query : getSearchStrings(true)) {
            Integer match0 = SequenceMatcher.fuzzyMatch(query, rec.getReadString().toUpperCase(), mismatchesAllowed);
            if (match0 != null) {
                int pos = rec.getReferencePositionAtReadPosition(match0);  //this is equal to one base prior to the 1-based start

                int i = 1;
                while (pos == 0 && (match0 - i) > 1) {
                    pos = rec.getReferencePositionAtReadPosition(match0 - i);
                    i++;
                }

                if (pos > 0) {
                    IntegrationSiteMapper.JunctionMatch jm = new IntegrationSiteMapper.JunctionMatch(log, id, this, rec.getContig(), pos, true, rec, maxRecordsToStore);
                    matches.put(jm.getKey(), jm);
                }
            }
        }

        return matches;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<String> getSearchStrings() {
        return searchStrings;
    }

    public void setSearchStrings(List<String> searchStrings) {
        this.searchStrings = searchStrings;
    }

    public List<String> getSearchStringsRC() {
        return searchStringsRC;
    }

    public void setSearchStringsRC(List<String> searchStringsRC) {
        this.searchStringsRC = searchStringsRC;
    }

    public int getMismatchesAllowed() {
        return mismatchesAllowed;
    }

    public boolean isInvertHitOrientation() {
        return invertHitOrientation;
    }

    public void setInvertHitOrientation(boolean invertHitOrientation) {
        this.invertHitOrientation = invertHitOrientation;
    }
}
