package com.github.discvrseq.util.help;

import jdk.javadoc.doclet.DocletEnvironment;
import org.broadinstitute.barclay.help.DocWorkUnit;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.help.GSONWorkUnit;
import org.broadinstitute.hellbender.utils.help.GATKDocWorkUnit;
import org.broadinstitute.hellbender.utils.help.GATKHelpDocWorkUnitHandler;
import org.broadinstitute.hellbender.utils.help.GATKHelpDoclet;

import javax.lang.model.element.Element;
import java.util.List;
import java.util.Map;

@SuppressWarnings("removal")
public class DISCVRSeqHelpDoclet extends GATKHelpDoclet {
    public DISCVRSeqHelpDoclet() {

    }

    public static boolean processDocs(final DocletEnvironment docletEnv) {
        return new DISCVRSeqHelpDoclet().run(docletEnv);
    }

    @Override
    public DocWorkUnit createWorkUnit(
            final Element classElement,
            final Class<?> clazz,
            final DocumentedFeature documentedFeature)
    {
        return new GATKDocWorkUnit(
                new GATKHelpDocWorkUnitHandler(this),
                classElement,
                clazz,
                documentedFeature);
    }

    @Override
    protected GSONWorkUnit createGSONWorkUnit(
            final DocWorkUnit workUnit,
            final List<Map<String, String>> groupMaps,
            final List<Map<String, String>> featureMaps)
    {
        return super.createGSONWorkUnit(workUnit, groupMaps, featureMaps);
    }
}
