package com.github.discvrseq.util;

import com.opencsv.*;
import htsjdk.samtools.util.IOUtil;

import java.io.File;
import java.nio.file.Path;

public class CsvUtils {
    public static ICSVWriter getTsvWriter(File output) {
        return getTsvWriter(output.toPath());
    }

    public static ICSVWriter getTsvWriter(Path output) {
        return new CSVWriterBuilder(IOUtil.openFileForBufferedUtf8Writing(output)).withSeparator('\t').withQuoteChar(CSVWriter.NO_QUOTE_CHARACTER).build();
    }

    public static CSVReader getTsvReader(File input) {
        return new CSVReaderBuilder(IOUtil.openFileForBufferedUtf8Reading(input)).withCSVParser(new CSVParserBuilder().withSeparator('\t').build()).build();
    }
}
