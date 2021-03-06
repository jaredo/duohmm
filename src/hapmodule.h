//$Id$

#pragma once
#include <vector>
#include <list>
#include <set>
#include <iterator>
#include <sstream>
#include <string>
#include <map>
#include "utilityfunc.h"
#include <cstring>

#ifdef SHAPEIT
#include <duohmm/utilityfunc.h>
#include <containers/genhap_set.h>
#include <output/haplotype_writer.h>

//S2 headers
#include <output/filter_writer.h>
#include <containers/genhap_set.h>
#endif

using namespace std;

class Haplotypes {
 public:
  Haplotypes(string filename);
  Haplotypes(const Haplotypes& h);
  ~Haplotypes();
  int writeHaps(string fname);
  string input_file;
  vector<bool> *getHap(string id);  
  vector<string> ids;
  vector<string> rsid1,rsid2,ref,alt;  
  vector<int> positions;
  map<string,int> idlook;//stores index of samples
  vector<float> cM;
  int nsnp,nhap,K,min_canopy_size,max_canopy_size;
  vector<bool> *H;  // haps raw/compressed
  vector<int> allsamples;
  bool isMissing(string & id,int index);
  vector<bool> *getMissing(string & id);
  vector < vector<bool> > _missing;
  int setMissing(const string & fname);
  vector<float> Lmatch,Rmatch;
  float impute(int hap_index,int pos_index);
  int imputeAll();
  int imputeSample(string const & id);
  void impute(string const & id,int pos_index);
      
#ifdef SHAPEIT
  Haplotypes(filter_writer & F, genhap_set & GH);
  int getSHAPEIT2(filter_writer & F, genhap_set & GH);
  vector < haplotype_index > orderI;
  filter_writer * _F;
  genhap_set * _GH;
#endif

};

