CXX ?= g++
CFLAGS ?= -O3

CINCLUDE=-I../include
HEADERS=baseheader.h multi-array.h util.h Multiseq.h nucleotide.h Aligner.h Consensus.h Progress.h

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS += -DDEBUG
else
	CFLAGS += -DNDEBUG
endif

TARGET: libcommon-bio.a

all: $(TARGET)

libcommon-bio.a: Consensus.o Aligner.o nucleotide.o Multiseq.o util.o baseheader.o Progress.o
	ar rcs $@ $^

Progress.o: Progress.cpp Progress.h
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

Consensus.o: Consensus.cpp Consensus.h nucleotide.h baseheader.h util.h
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

Aligner.o: Aligner.cpp Aligner.h util.h multi-array.h baseheader.h
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

nucleotide.o: nucleotide.cpp nucleotide.h baseheader.h
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

Multiseq.o: Multiseq.cpp Multiseq.h util.h baseheader.h
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

util.o: util.cpp util.h baseheader.h
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

baseheader.o: baseheader.cpp baseheader.h
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o
	
