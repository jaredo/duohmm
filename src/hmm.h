#pragma once
#include "hapmodule.h"
#include "utilityfunc.h"

class geneticMap{
 public:
  geneticMap(string fname);
  int interpolate(vector<int> & positions,vector<double> & output);
  int nsnp;
  double interpolate(int position);
 private:
  vector<int> pos;
  vector<double> cM;// cM[t] is the genetic distance between snp t and t+1
};

class DuoHMM{
 public:
  DuoHMM(vector<int> & positions, geneticMap & gm);
  int nsnp,K;

  double error,switch1,switch2;//parameters
  unsigned char **parent,**child;//pointers to duo haps

  //recombination maps
  double male_multiplier,female_multiplier,genetic_length;
  vector <double> male_rho;
  vector <double> female_rho;
  vector <double> male_norho;
  vector <double> female_norho;
  vector<double> cM;
  double *rho;

  int setHaps(unsigned char **parental_haplotypes,unsigned char **child_haplotypes,string sex);
  int EM(int niteration);

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


