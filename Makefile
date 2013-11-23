CC = g++
BOOSTFLAGS = -lboost_iostreams -lboost_random -lboost_program_options
#CFLAGS = -O3  -std=gnu++0x $(BOOSTFLAGS)
CFLAGS = -O3  -Wall $(BOOSTFLAGS)
#CFLAGS = -g -Wall  $(BOOSTFLAGS)

duohmm: pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o  src/duohmm.cpp 
	$(CC)  src/duohmm.cpp pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o -o duohmm $(CFLAGS)
test: pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o src/test.cpp
	$(CC)  src/test.cpp pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o -o test $(CFLAGS)
pedhap.o: src/pedhap.h src/pedhap.cpp utilityfunc.o pedigree.o hmm_duo.o hmm_trio.o
	$(CC) -c src/pedhap.cpp $(CFLAGS) 
hapmodule.o: src/hapmodule.cpp src/hapmodule.h utilityfunc.o
	$(CC) -c src/hapmodule.cpp $(CFLAGS) 
utilityfunc.o:  src/utilityfunc.cpp src/utilityfunc.h
	$(CC) -c src/utilityfunc.cpp $(CFLAGS)
pedigree.o:  src/pedigree.h src/pedigree.cpp src/utilityfunc.cpp src/utilityfunc.h
	$(CC) -c src/pedigree.cpp $(CFLAGS)
hmm_duo.o:  src/hmm.h src/hmm_duo.cpp 
	$(CC) -c src/hmm_duo.cpp $(CFLAGS)
hmm_trio.o:  src/hmm.h src/hmm_trio.cpp
	$(CC) -c src/hmm_trio.cpp $(CFLAGS)

all: duohmm test

clean:
	rm *.o
	rm duohmm test
