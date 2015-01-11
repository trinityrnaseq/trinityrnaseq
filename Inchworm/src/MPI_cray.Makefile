SOURCES=Fasta_entry.cpp sequenceUtil.cpp \
                   MPIinchworm.cpp KmerCounter.cpp string_util.cpp \
                   Fasta_reader.cpp stacktrace.cpp argProcessor.cpp


# module swap PrgEnv-cray PrgEnv-intel

CC=CC
CFLAGS=-openmp -Wno-deprecated -O3
LDFLAGS=-openmp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MPIinchworm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@
