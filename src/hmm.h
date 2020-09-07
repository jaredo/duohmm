//$Id$

#pragma once
#include "hapmodule.h"
#include "utilityfunc.h"

#define DEBUG 0

class geneticMap{
 public:
#ifdef SHAPEIT
  geneticMap(genhap_set & GH);
#endif
  geneticMap(string fname);
  geneticMap();
  int interpolate(vector<int> & positions,vector<float> & output);
  int nsnp;
  float interpolate(int position);
 private:
  vector<float> pos;
  vector<float> cM;// cM[t] is the genetic distance between snp t and t+1
};

class DuoHMM{
 public:
  DuoHMM(vector<int> & positions, geneticMap & gm);
  int nsnp,K,NITERATION;

  float error,switch1,switch2;//parameters
  vector<bool> *parent,*child;//pointers to duo haps

  //recombination maps
  float male_multiplier,female_multiplier,male_length,female_length,genetic_length,multi;
  vector <float> male_rho;
  vector <float> female_rho;
  vector <float> male_norho;
  vector <float> female_norho;
  vector <float> cM;
  vector <float> recombinationMap;
  vector<float> genError;

  float *rho;

  int setHaps(vector<bool> *parental_haplotypes,vector<bool> *child_haplotypes,string sex);
  void setIterations(int n);

  int EM(int niteration);
  int estimateRecombination();


  //F-B VARIABLES
  vector < vector<float> > alpha;
  vector < vector<float> > beta;
  vector <float> scale;
  int forward();
  int backward();
  vector < vector<float> > posterior;
  vector < vector< vector<float> > > trans_posterior;

  int estep();
  int mstep();

  //VITERBI VARIABLES
  vector< vector<unsigned char> > backtrack;
  int viterbi();
  vector<unsigned char> stateseq;
};


class TrioHMM {
 public:
  TrioHMM(vector<int> & positions, geneticMap & gm);
  int nsnp,K,NITERATION;

  float error,switch_child,switch_mum,switch_dad;//parameters
  vector<bool> *dad, *mum, *child;//pointers to duo haps

  //recombination maps
  float male_multiplier,female_multiplier,male_length,female_length;
  vector <float> male_rho;
  vector <float> female_rho;
  vector <float> male_norho;
  vector <float> female_norho;
  vector<float> cM;
  vector <float> recombinationMat;
  vector <float> recombinationPat;
  int setHaps(vector<bool> *dadptr,vector<bool> *mumptr,vector<bool> *childptr);
  int EM(int niteration);
  void setIterations(int n);

  //F-B VARIABLES
  vector < vector<float> > alpha;
  vector < vector<float> > beta;
  vector <float> scale;
  int forward();
  int backward();
  int estimateRecombination();
  int estimateRecombinationDad();
  int estimateRecombinationMum();
  vector < vector<float> > posterior;
  vector < vector< vector<float> > > trans_posterior;

  int estep();
  int mstep();
  vector<float> genError;

  //VITERBI VARIABLES
  vector< vector<unsigned char> > backtrack;
  int viterbi();
  vector<unsigned char> stateseq;
};


