#include "utilityfunc.h"


int argmax(vector<double> & x) {
  int maxind=0;
  for(int i=1;i<x.size();i++) 
    if(x[i]>x[maxind]) maxind=i;
  return(maxind);
}


unsigned char**uc2hap() {
  int n = 256;
  unsigned char **ret = newMatrix<unsigned char>(n,8);
  for(unsigned int i=0;i<n;i++) {
    bitset<8> bs(i);
    for(int j=0;j<8;j++) ret[i][j] = bs[j];
  }
  return(ret);

}

unsigned char** buildLookup() {
  int n = 256;
  unsigned char **ret = newMatrix<unsigned char>(n,n);
  for(unsigned int i=0;i<n;i++) {
    bitset<8> bs1(i);
    for(unsigned int j=0;j<n;j++) {
      bitset<8> bs2(j);
      ret[i][j] = (unsigned char)(bs1^bs2).count();
    }
  }
  return(ret);
}

int which_max(int *x,int n) {
  int maxidx = 0;
  int maxval = x[maxidx];
  for(int i=0;i<n;i++) {
    if(x[i]>maxval) {
      maxidx = i;
      maxval = x[i];
    }
  }
  return maxidx;
}

bool fileexists(string fname){
  ifstream ifile(fname.c_str());
  return ifile.is_open();
}
