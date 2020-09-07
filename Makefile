CXX = g++
$(shell mkdir -p bin/ obj/) 

ifdef BOOST_ROOT
	LFLAGS = -lz $(BOOST_ROOT)/lib/libboost_iostreams.a $(BOOST_ROOT)/lib/libboost_program_options.a 
	ifneq (,$(wildcard $(BOOST_ROOT)/lib64/libboost_iostreams.a))
		LFLAGS = -lz $(BOOST_ROOT)/lib64/libboost_iostreams.a $(BOOST_ROOT)/lib64/libboost_program_options.a
	endif
	IFLAGS = -I$(BOOST_ROOT)/include
else
	LFLAGS = -lz -lboost_iostreams -lboost_program_options
endif

CFLAGS = -O3 $(IFLAGS) -Wall -std=c++11 

VERSION := $(shell git describe  --always --tags)

all: bin/duohmm
debug: CFLAGS = -g -Wall $(IFLAGS) -std=c++11 
debug: all

OBJS=obj/pedhap.o obj/hapmodule.o obj/utilityfunc.o obj/pedigree.o obj/hmm_duo.o obj/hmm_trio.o  
obj/%.o: src/%.cpp
	$(CXX) -o $@ -c $< $(CFLAGS) $(IFLAGS)
bin/duohmm: $(OBJS) src/duohmm.cpp 
	$(CXX) -DVERSION=\"$(VERSION)\" $(OBJS) src/duohmm.cpp -o $@ $(CFLAGS)  $(LFLAGS) -Ilib/
bin/sparseimpute: $(OBJS) src/sparseimpute.cpp 
	$(CXX) -DVERSION=\"$(VERSION)\" $(OBJS) src/sparseimpute.cpp -o $@ $(CFLAGS)  $(LFLAGS) -Ilib/	
bin/unit_test: $(OBJS) src/test.cpp
	$(CXX)  test.cpp $(OBJS) -o test $(CFLAGS)

test: bin/duohmm
	cd test/;bash -e test.sh
clean:
	rm -rf obj/ bin/ test/observed.*
