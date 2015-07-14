ifndef BOOST_ROOT
    $(error BOOST_ROOT is undefined)
endif

CC = g++
#LFLAGS = -L$(BOOST_ROOT)/lib -lz -lboost_iostreams  -lboost_program_options
LFLAGS = -lz $(BOOST_ROOT)/lib/libboost_iostreams.a $(BOOST_ROOT)/lib/libboost_program_options.a
IFLAGS =  -I$(BOOST_ROOT)/include

VERSION := $(shell git describe --dirty --always --tags)

release: CFLAGS = -O3 $(IFLAGS)
release: duohmm

debug: CFLAGS = -g -Wall $(IFLAGS)
debug: duohmm 

duohmm: pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o  duohmm.cpp 
	$(CC) -DVERSION=\"$(VERSION)\" duohmm.cpp pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o -o duohmm $(CFLAGS)  $(LFLAGS)
test: pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o test.cpp
	$(CC)  test.cpp pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o -o test $(CFLAGS)
pedhap.o: pedhap.h pedhap.cpp utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o
	$(CC) -c pedhap.cpp $(CFLAGS) 
hapmodule.o: hapmodule.cpp hapmodule.h utilityfunc.o
	$(CC) -c hapmodule.cpp $(CFLAGS) 
utilityfunc.o:  utilityfunc.cpp utilityfunc.h
	$(CC) -c utilityfunc.cpp $(CFLAGS)
pedigree.o:  pedigree.h pedigree.cpp utilityfunc.cpp utilityfunc.h
	$(CC) -c pedigree.cpp $(CFLAGS)
hmm_duo.o:  hmm.h hmm_duo.cpp 
	$(CC) -c hmm_duo.cpp $(CFLAGS)
hmm_trio.o:  hmm.h hmm_trio.cpp
	$(CC) -c hmm_trio.cpp $(CFLAGS)

all: duohmm test

clean:
	rm *.o
	rm duohmm test
