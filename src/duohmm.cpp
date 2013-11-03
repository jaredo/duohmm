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
#include <boost/progress.hpp>

using namespace std;

int main(int argc,char **argv) {

	if(argc!=3) {
		cout << "Usage: duohmm hapfile famfile --corrected output --recombination output --genotypingError output" << endl;
		return(0);
	}

	cout << "Reading haplotypes from " << argv[1] <<".haps" << endl;
	Haplotypes haps(argv[1]);

	cout << "Reading pedigree information from " << argv[2] << endl;

	pedigree p(argv[2]);
	return 0;
}
