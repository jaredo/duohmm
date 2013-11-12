CC = g++
BOOSTFLAGS = -lboost_iostreams -lboost_random -lboost_program_options
#CFLAGS = -O3  -std=gnu++0x $(BOOSTFLAGS)
CFLAGS = -O3  $(BOOSTFLAGS)
#CFLAGS = -g -Wall  $(BOOSTFLAGS)

test: pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm.o src/test.cpp
	$(CC)  src/test.cpp pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm.o -o test $(CFLAGS)
duohmm: pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm.o src/duohmm.cpp 
	$(CC)  src/duohmm.cpp pedhap.o hapmodule.o utilityfunc.o pedigree.o hmm.o -o duohmm $(CFLAGS)
pedhap.o: src/pedhap.h src/pedhap.cpp utilityfunc.o pedigree.o hmm.o
	$(CC) -c src/pedhap.cpp $(CFLAGS) 
hapmodule.o: src/hapmodule.cpp src/hapmodule.h utilityfunc.o
	$(CC) -c src/hapmodule.cpp $(CFLAGS) 
utilityfunc.o:  src/utilityfunc.cpp src/utilityfunc.h
	$(CC) -c src/utilityfunc.cpp $(CFLAGS)
pedigree.o:  src/pedigree.h src/pedigree.cpp src/utilityfunc.cpp src/utilityfunc.h
	$(CC) -c src/pedigree.cpp $(CFLAGS)
hmm.o:  src/hmm.h src/hmm.cpp
	$(CC) -c src/hmm.cpp $(CFLAGS)

all: duohmm test

clean:
	rm *.o
	rm duohmm test
