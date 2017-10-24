#pragma once
#include "utilityfunc.h"
#include "pedigree.h"
#include "hmm.h"

#ifndef SHAPEIT
#define VERBOSE 1
#endif

#ifdef SHAPEIT
#include <globals.h>
#define VERBOSE 0
#endif

class pedhap {
 public:
  pedhap(string hap_filename,string pedigree_filename,string gm_filename);
  pedhap(string hap_filename,string gm_filename,int niteration=100);
  ~pedhap();
#ifdef SHAPEIT
  pedhap(filter_writer & F, genhap_set & GH,string header1,string header2,int niteration=100);
#endif
  //  pedhap(string hap_filename,string pedigree_filename);


  int phase(string child);
  int correct();
  int minRecombinant(string parent);
  int genotypingError(string outfile);
  int recombinationMap(string outfile);

  int nsnp,NITERATION;
  map< pair<string,string>, vector<unsigned char> > stateseq;

  geneticMap *gm;
  DuoHMM *duo;
  TrioHMM *trio;
  Haplotypes *haps;
  pedigree *ped;
};
