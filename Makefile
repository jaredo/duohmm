CC = g++
BOOSTFLAGS = -lboost_iostreams -lboost_random -lboost_program_options
CFLAGS = -O3  -std=gnu++0x $(BOOSTFLAGS)

kfinder: hapmodule.o utilityfunc.o pedigree.o src/duohmm.cpp
	$(CC)  src/duohmm.cpp hapmodule.o utilityfunc.o pedigree.o -o duohmm $(CFLAGS)
hapmodule.o: src/hapmodule.cpp src/hapmodule.h utilityfunc.o
	$(CC) -c src/hapmodule.cpp $(CFLAGS) 
utilityfunc.o:  src/utilityfunc.cpp src/utilityfunc.h
	$(CC) -c src/utilityfunc.cpp $(CFLAGS)
pedigree.o:  src/pedigree.h src/pedigree.cpp src/utilityfunc.cpp src/utilityfunc.h
	$(CC) -c src/pedigree.cpp $(CFLAGS)

all: duohmm
clean:
	rm *.o
	rm duohmm
