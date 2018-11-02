CC ?= gcc
CXX ?= g++
CFLAGS ?= -O3
CINCLUDE=-I../include
LIBDIR=../lib

TARGET=findmaxmat.o prenuc.o mgaps.o postnuc.o delta-filter.o nucmer.o

all: $(TARGET)

nucmer.o: nucmer.cpp findmaxmat.o mgaps.o postnuc.o delta-filter.o
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@
	ar rcs libmummer.a nucmer.o findmaxmat.o mgaps.o postnuc.o delta-filter.o tigrinc.o sw_align.o delta.o

findmaxmat.o: findmaxmat.c
	$(CC) $(CFLAGS) -c $< $(CINCLUDE) -o $@

prenuc.o: prenuc.cpp
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@

mgaps.o: mgaps.cpp tigrinc.o
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@

postnuc.o: postnuc.cpp tigrinc.o sw_align.o
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@

delta-filter.o: delta-filter.cpp delta.o tigrinc.o
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@

tigrinc.o: tigrinc.cpp
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@

sw_align.o: sw_align.cpp
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@

delta.o: delta.cpp
	$(CXX) $(CFLAGS) -c $< $(CINCLUDE) -o $@


clean:
	rm -f *.o


