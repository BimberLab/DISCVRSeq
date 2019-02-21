package com.github.discvrseq.util.help;

import com.sun.javadoc.ClassDoc;
import com.sun.javadoc.RootDoc;
import org.broadinstitute.barclay.help.DocWorkUnit;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.help.GSONWorkUnit;
import org.broadinstitute.hellbender.utils.help.GATKDocWorkUnit;
import org.broadinstitute.hellbender.utils.help.GATKHelpDocWorkUnitHandler;
import org.broadinstitute.hellbender.utils.help.GATKHelpDoclet;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public class DISCVRSeqHelpDoclet extends GATKHelpDoclet {
    public DISCVRSeqHelpDoclet() {

    }

    /**
     * Create a doclet of the appropriate type and generate the FreeMarker templates properties.
     * @param rootDoc
     * @throws IOException
     */
    public static boolean start(final RootDoc rootDoc) throws IOException {
        return new DISCVRSeqHelpDoclet().startProcessDocs(rootDoc);
    }

    /**
     * @return Create and return a DocWorkUnit-derived object to handle documentation
     * for the target feature(s) represented by documentedFeature.
     *
     * @param documentedFeature DocumentedFeature annotation for the target feature
     * @param classDoc javadoc classDoc for the target feature
     * @param clazz class of the target feature
     * @return DocWorkUnit to be used for this feature
     */
    @Override
    protected DocWorkUnit createWorkUnit(
            final DocumentedFeature documentedFeature,
            final ClassDoc classDoc,
            final Class<?> clazz)
    {
        return new GATKDocWorkUnit(
                new GATKHelpDocWorkUnitHandler(this),
                documentedFeature,
                classDoc,
                clazz);
    }

    /**
     * Create a GSONWorkUnit-derived object that holds our custom data. This method should create the object, and
     * propagate any custom javadoc tags from the template map to the newly created GSON object; specifically
     * "walkertype", which is pulled from a custom javadoc tag.
     *
     * @param workUnit work unit for which a GSON object is required
     * @param groupMaps
     * @param featureMaps
     * @return a GSONWorkUnit-derived object for this work unit, populated with any custom values
     */
    @Override
    protected GSONWorkUnit createGSONWorkUnit(
            final DocWorkUnit workUnit,
            final List<Map<String, String>> groupMaps,
            final List<Map<String, String>> featureMaps)
    {
        return super.createGSONWorkUnit(workUnit, groupMaps, featureMaps);
    }
}
