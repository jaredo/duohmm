//============================================================================
// Name        : duohmm.cpp
// Author      : Jared O'Connell
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "hapmodule.h"
#include "pedigree.h"
#include "hmm.h"
#include "pedhap.h"
#include <boost/progress.hpp>

using namespace std;

int main(int argc,char **argv) {

  if(argc!=4) {
    cout << "Usage: duohmm hapfile famfile genetic_map --corrected output --recombination output --genotypingError output" << endl;
    return(0);
  }

  pedhap ph(argv[1],argv[2],argv[3]);
  ph.correct();
  ph.haps->writeHaps(argv[1]);
  return 0;
  
}
