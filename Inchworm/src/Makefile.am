AUTOMAKE_OPTIONS = foreign

bin_PROGRAMS = inchworm cigar_tweaker pull_reads_with_kmers FastaToDeBruijn fastaToKmerCoverageStats

inchworm_SOURCES = Fasta_entry.cpp IRKE_run.cpp  sequenceUtil.cpp \
                   IRKE.cpp KmerCounter.cpp string_util.cpp \
                   Fasta_reader.cpp IRKE.hpp stacktrace.cpp argProcessor.cpp


cigar_tweaker_SOURCES = cigar_tweaker.cpp SAM_entry.cpp sequenceUtil.cpp SAM_reader.cpp Cigar.cpp string_util.cpp stacktrace.cpp Fasta_reader.cpp Fasta_entry.cpp


pull_reads_with_kmers_SOURCES = pull_reads_with_kmers.cpp Fasta_reader.cpp Fasta_entry.cpp sequenceUtil.cpp KmerCounter.cpp string_util.cpp stacktrace.cpp argProcessor.cpp

FastaToDeBruijn_SOURCES = FastaToDeBruijn.cpp argProcessor.cpp Fasta_reader.cpp Fasta_entry.cpp sequenceUtil.cpp string_util.cpp stacktrace.cpp DeBruijnGraph.cpp

fastaToKmerCoverageStats_SOURCES = fastaToKmerCoverageStats.cpp argProcessor.cpp Fasta_reader.cpp Fasta_entry.cpp sequenceUtil.cpp string_util.cpp stacktrace.cpp KmerCounter.cpp 
