//$Id$

#pragma once
#include <utility> 
#include <vector>
#include <bitset>
#include <assert.h>

#include <cmath> 
#include <sys/stat.h>
#include <fstream>
#include <iostream>

using namespace std;

unsigned char**uc2hap() ;

unsigned char** buildLookup();

int which_max(int *x,int n);

template <class T>
bool comparator ( const pair<T,int>& l, const pair<T,int>& r)  { return l.first < r.first; };

template <class T>
vector<int> argsort(vector<T> *list) {

  vector<pair<T,int> > tosort(list->size());
  for(int i=0;i<list->size();i++) {
    tosort[i] = make_pair( (*list)[i] ,i );
  }
  sort(tosort.begin(),tosort.end());//,comparator);
  vector<int> ret(tosort.size());
  for(int i=0;i<tosort.size();i++) {
    ret[i] = tosort[i].second;
    //    cout << tosort[i].first << "\t" << tosort[i].second << endl;
  }
  return(ret);
}

inline int ham(unsigned char*x1,unsigned char*x2,int n,unsigned char**look) {
  int d = 0;
  for(int i=0;i<n;i++)
    d += look[x1[i]][x2[i]];
  return(d);
}

template <class T>
T **newMatrix(int nrow,int ncol) {
  T **ret = new T*[nrow];
  for(int i=0;i<nrow;i++)
    ret[i] = new T[ncol];  
  return(ret);
}


template <class T>
void delMatrix(T **mat,int nrow) {
  for(int i=0;i<nrow;i++)
    delete[] mat[i];
  delete[] mat;
}

template <class T>
void printMatrix(T **H,int nrow,int ncol) {
  for(int i=0;i<nrow;i++) {
    for(int j=0;j<ncol;j++)
      cout << (unsigned int) H[i][j] << "\t";
    cout << endl;
  }
}

bool fileexists(string fname);

int argmax(vector<float> & x);

bool is_mendel_consistent(int k,int d,int m);
