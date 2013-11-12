#pragma once
#include "utilityfunc.h"
#include "pedigree.h"
#include "hmm.h"


class pedhap {
 public:
  pedhap(string hap_filename,string pedigree_filename,string gm_filename);
  int phase(string parent,string child);
  int correct();
  int minRecombinant(string parent);
  //  int genotypingError();
  //  int recombinationMap();
  int nsnp;
  map< pair<string,string>, vector<unsigned char> > stateseq;

  geneticMap *gm;
  DuoHMM *duo;
  Haplotypes *haps;
  pedigree *ped;
};
