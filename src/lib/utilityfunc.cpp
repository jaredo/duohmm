#include "utilityfunc.h"


int argmax(vector<double> & x) {
  int maxind=0;
  for(size_t i=1;i<x.size();i++) 
    if(x[i]>x[maxind]) maxind=i;
  return(maxind);
}


unsigned char**uc2hap() {
  size_t n = 256;
  unsigned char **ret = newMatrix<unsigned char>(n,8);
  for(size_t i=0;i<n;i++) {
    bitset<8> bs(i);
    for(int j=0;j<8;j++) ret[i][j] = bs[j];
  }
  return(ret);

}

unsigned char** buildLookup() {
  size_t n = 256;
  unsigned char **ret = newMatrix<unsigned char>(n,n);
  for(size_t i=0;i<n;i++) {
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

bool is_mendel_consistent(int k,int d,int m)
{
    if(m==0&&d==0)
    {
	if(k==0) return(true);
    }
    if((m==0&&d==1)||(m==1&&d==0))
	if(k!=2) return(true);
    if((m==0&&d==2)||(m==2&&d==0))
	if(k==1) return(true);
    if(m==1&&d==1)
	return(true);
    if((m==1&&d==2)||(m==2&&d==1))
	if(k!=0) return(true);      
    if(m==2&&d==2)
	if(k==2) return(true);

    return(false);
}
