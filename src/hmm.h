#pragma once
#include "hapmodule.h"
#include "utilityfunc.h"

#define NITERATION 100
#define DEBUG 0


class geneticMap{
 public:
  geneticMap(string fname);
  geneticMap();
  int interpolate(vector<int> & positions,vector<double> & output);
  int nsnp;
  double interpolate(int position);
 private:
  vector<double> pos;
  vector<double> cM;// cM[t] is the genetic distance between snp t and t+1
};

class DuoHMM{
 public:
  DuoHMM(vector<int> & positions, geneticMap & gm);
  int nsnp,K;

  double error,switch1,switch2;//parameters
  unsigned char **parent,**child;//pointers to duo haps

  //recombination maps
  double male_multiplier,female_multiplier,male_length,female_length,genetic_length,multi;
  vector <double> male_rho;
  vector <double> female_rho;
  vector <double> male_norho;
  vector <double> female_norho;
  vector <double> cM;
  vector <double> recombinationMap;
  vector<double> genError;

  double *rho;

  int setHaps(unsigned char **parental_haplotypes,unsigned char **child_haplotypes,string sex);
  int EM(int niteration);
  int estimateRecombination();


  //F-B VARIABLES
  vector < vector<double> > alpha;
  vector < vector<double> > beta;
  vector <double> scale;
  int forward();
  int backward();
  vector < vector<double> > posterior;
  vector < vector< vector<double> > > trans_posterior;

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
  int nsnp,K;

  double error,switch_child,switch_mum,switch_dad;//parameters
  unsigned char **dad, **mum, **child;//pointers to duo haps

  //recombination maps
  double male_multiplier,female_multiplier,male_length,female_length;
  vector <double> male_rho;
  vector <double> female_rho;
  vector <double> male_norho;
  vector <double> female_norho;
  vector<double> cM;
  vector <double> recombinationMat;
  vector <double> recombinationPat;
  int setHaps(unsigned char **dadptr,unsigned char **mumptr,unsigned char **childptr);
  int EM(int niteration);

  //F-B VARIABLES
  vector < vector<double> > alpha;
  vector < vector<double> > beta;
  vector <double> scale;
  int forward();
  int backward();
  int estimateRecombination();
  vector < vector<double> > posterior;
  vector < vector< vector<double> > > trans_posterior;

  int estep();
  int mstep();
  vector<double> genError;

  //VITERBI VARIABLES
  vector< vector<unsigned char> > backtrack;
  int viterbi();
  vector<unsigned char> stateseq;
};


