package com.github.discvrseq.util;

import htsjdk.samtools.*;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

public class SamUtils {
    public static File ensureQuerySorted(File bam, Path genomeFasta, Logger logger) {
        SamReaderFactory fact = SamReaderFactory.makeDefault();
        fact.validationStringency(ValidationStringency.SILENT);
        fact.referenceSequence(genomeFasta);

        SAMFileHeader.SortOrder so = fact.getFileHeader(bam).getSortOrder();
        if (so == SAMFileHeader.SortOrder.queryname) {
            logger.info("BAM is already query sorted, no need to sort");
            return bam;
        }

        logger.info("Sorting BAM in queryName order");
        try (SamReader reader = SamReaderFactory.makeDefault().referenceSequence(genomeFasta).open(bam)) {
            reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.queryname);

            File querySorted = File.createTempFile(bam.getName(), ".querySorted.bam").toPath().normalize().toFile();
            logger.info("writing to file: " + querySorted.getPath());

            try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, querySorted)) {
                for (final SAMRecord rec : reader) {
                    writer.addAlignment(rec);
                }
            }

            return querySorted;
        }
        catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }
    }
}
