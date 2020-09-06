ifndef BOOST_ROOT
    $(error BOOST_ROOT is undefined)
endif

CXX = g++
LFLAGS = -lz -lboost_iostreams -lboost_program_options -L$(BOOST_ROOT)/lib 
LFLAGS = -lz $(BOOST_ROOT)/lib/libboost_iostreams.a $(BOOST_ROOT)/lib/libboost_program_options.a
IFLAGS = -I$(BOOST_ROOT)/include
CFLAGS = -O3 $(IFLAGS) -Wall -std=c++11 

VERSION := asdf#$(shell git describe  --always --tags)

all: duohmm
debug: CFLAGS = -g -Wall $(IFLAGS) -std=c++11 
debug: all

OBJS=obj/pedhap.o obj/hapmodule.o obj/utilityfunc.o obj/pedigree.o obj/hmm_duo.o obj/hmm_trio.o  
obj/:
	mkdir obj
obj/%.o: lib/%.cpp obj/
	$(CXX) -o $@ -c $< $(CFLAGS) $(IFLAGS)
duohmm: $(OBJS) duohmm.cpp 
	$(CXX) -DVERSION=\"$(VERSION)\" $(OBJS) duohmm.cpp -o $@ $(CFLAGS)  $(LFLAGS) -Ilib/
sparseimpute: $(OBJS) sparseimpute.cpp 
	$(CXX) -DVERSION=\"$(VERSION)\" $(OBJS) sparseimpute.cpp -o $@ $(CFLAGS)  $(LFLAGS) -Ilib/	
unit_test: $(OBJS) test.cpp
	$(CXX)  test.cpp $(OBJS) -o test $(CFLAGS)
obj/%.o: lib/utilityfunc.cpp lib/utilityfunc.h

test: duohmm
	cd ../test/;bash -e test.sh
clean:
	rm -rf obj/
	rm -f duohmm unit_test ../test/observed.*
