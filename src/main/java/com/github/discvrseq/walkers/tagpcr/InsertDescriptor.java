package com.github.discvrseq.walkers.tagpcr;

import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class InsertDescriptor {
    private String name;
    private SequenceDescriptor insertUpstreamRegion;
    private SequenceDescriptor insertDownstreamRegion;

    private List<SequenceDescriptor> internalPrimers;
    private List<String> backboneSearchStrings;
    private List<String> backboneSearchStringsRC = null;
    private int backboneSearchEditDistance = 2;

    private List<InsertJunctionDescriptor> junctions;

    public InsertDescriptor() {

    }

    public String getFeatureLabel(TagPcrSummary.END_TYPE type) {
        switch (type) {
            case upstream:
                return insertUpstreamRegion.getName();
            case downstream:
                return insertDownstreamRegion.getName();
            default:
                throw new IllegalArgumentException("Unknown type");
        }
    }

    public String getName() {
        return name;
    }

    public SequenceDescriptor getInsertUpstreamRegion() {
        return insertUpstreamRegion;
    }

    public SequenceDescriptor getInsertDownstreamRegion() {
        return insertDownstreamRegion;
    }

    public List<SequenceDescriptor> getInternalPrimers() {
        return internalPrimers;
    }

    public List<InsertJunctionDescriptor> getJunctions() {
        return junctions;
    }

    public int getBackboneSearchEditDistance() {
        return backboneSearchEditDistance;
    }

    public void setBackboneSearchEditDistance(int backboneSearchEditDistance) {
        this.backboneSearchEditDistance = backboneSearchEditDistance;
    }

    public List<String> getAllBackboneSearchStrings() {
        if (backboneSearchStrings == null) {
            return Collections.emptyList();
        }

        List<String> ret = new ArrayList<>(backboneSearchStrings);
        if (backboneSearchStringsRC == null) {
            backboneSearchStringsRC = new ArrayList<>();
            backboneSearchStrings.forEach(x -> {
                backboneSearchStringsRC.add(SequenceUtil.reverseComplement(x));
            });
        }

        ret.addAll(backboneSearchStringsRC);

        return ret;
    }

    public List<String> getBackboneSearchStrings() {
        return backboneSearchStrings;
    }

    public void setBackboneSearchStrings(List<String> backboneSearchStrings) {
        this.backboneSearchStrings = backboneSearchStrings;
    }
}
