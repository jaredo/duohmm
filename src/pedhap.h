#pragma once
#include "utilityfunc.h"
#include "pedigree.h"
#include "hmm.h"


class pedhap {
 public:
  pedhap(string hap_filename,string pedigree_filename,string gm_filename);
  pedhap(string hap_filename,string gm_filename);
  //  pedhap(string hap_filename,string pedigree_filename);


  int phase(string child);
  int correct();
  int minRecombinant(string parent);
  int genotypingError(string outfile);
  int recombinationMap(string outfile);

  int nsnp;
  map< pair<string,string>, vector<unsigned char> > stateseq;

  geneticMap *gm;
  DuoHMM *duo;
  TrioHMM *trio;
  Haplotypes *haps;
  pedigree *ped;
};
