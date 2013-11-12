#pragma once
#include <map>
#include <vector>
#include <set>
#include "utilityfunc.h"

struct individual{
  string fid;
  string id;
  string dad;
  string mum;
  string sex;
  int idx;
  int fidx;//pedigree index
  set<string> kids;
};

int buildFam(string id, map<string,individual> & sampleinfo, set<string> & fam);

class pedigree {
 public:
  pedigree(string fname,vector<string> & ids);
  int orderSamples(set<string> & ids,vector<string> & ordered_ids);
  map<string,individual> sampleinfo;
  vector<set<string> > pedigrees;
  int buildPeds();
};
