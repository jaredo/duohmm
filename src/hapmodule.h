#pragma once
#include <vector>
#include <list>
#include <set>
#include <unordered_set>
#include <iterator>
#include "utilityfunc.h"

#define DEBUG 0
using namespace std;

class Haplotypes {
 public:
  Haplotypes(string filename);
  Haplotypes(const Haplotypes& h);
  ~Haplotypes();

  int nsnp,nhap,K,min_canopy_size,max_canopy_size;
  unsigned char **H;  // haps raw/compressed
  vector<int> allsamples;
};
