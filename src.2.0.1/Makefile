CXX ?= g++
CC ?= gcc
MAKE ?= make
CFLAGS=-O3
CINCLUDE=-I./include
CLIB=-lemon -lpthread -L./lib -lcommon-bio -lmummer -lstree -lbase

COMMONBIO_HEADERS=./common-bioinfo/baseheader.h ./common-bioinfo/multi-array.h ./common-bioinfo/util.h ./common-bioinfo/Multiseq.h ./common-bioinfo/nucleotide.h ./common-bioinfo/Aligner.h ./common-bioinfo/Consensus.h ./common-bioinfo/Progress.h

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS += -DDEBUG
else
	CFLAGS += -DNDEBUG
endif

export CXX CC CFLAGS DEBUG

TARGET=smsc blasrm5_filter
OBJS=mummer_api.o threadpool.o path_cover.o
LIBS=./lib/libcommon-bio.a ./lib/libmummer.a ./lib/libstree.a ./lib/libbase.a

all: $(TARGET)

smsc: smsc.2.0.1.cpp $(OBJS) run_nucmer.sh
	$(CXX) $(CFLAGS) $(CINCLUDE) $< $(OBJS) -o $@ $(CLIB)

blasrm5_filter: blasrm5_filter.cpp
	$(CXX) $(CFLAGS) $^ -o $@

path_cover.o: path_cover.cpp path_cover.h
	$(CXX) $(CFLAGS) -c $< -o $@

mummer_api.o: mummer_api.cpp mummer_api.h
	mkdir -p lib
	cd libbasedir; $(MAKE) libbase.a
	cp libbasedir/libbase.a lib
	cd streesrc; $(MAKE) libstree.a
	cp streesrc/libstree.a lib
	cd mummer-lib; $(MAKE)
	cp mummer-lib/libmummer.a lib
	cd common-bioinfo; $(MAKE)
	cp common-bioinfo/libcommon-bio.a lib
	mkdir -p include/common-bio
	cp $(COMMONBIO_HEADERS) include/common-bio/
	$(CXX) $(CFLAGS) $(CINCLUDE) -c $< -o $@

threadpool.o: threadpool.c threadpool.h
	$(CC) -O3 -c -fPIC -D_REENTRANT -Wall -pedantic -I src -o $@ $<

clean:
	cd libbasedir; $(MAKE) clean
	cd streesrc; $(MAKE) clean
	cd mummer-lib; $(MAKE) clean
	cd common-bioinfo; $(MAKE) clean
	rm -rf include/common-bio
	rm -f lib/*.a
	rm -f *.o $(TARGET)
