#pragma once
#include <vector>
#include <list>
#include <set>
#include <iterator>
#include <sstream>
#include <string>
#include <map>
#include "utilityfunc.h"


using namespace std;

class Haplotypes {
 public:
  Haplotypes(string filename);
  Haplotypes(const Haplotypes& h);
  ~Haplotypes();
  int writeHaps(string fname);

  unsigned char **getHap(string id);  
  vector<string> ids;
  vector<string> rsid1,rsid2,ref,alt;  
  vector<int> positions;
  map<string,int> idlook;//stores index of samples
  vector<double> cM;
  int nsnp,nhap,K,min_canopy_size,max_canopy_size;
  unsigned char **H;  // haps raw/compressed
  vector<int> allsamples;

};
