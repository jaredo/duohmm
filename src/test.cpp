//============================================================================
// Name        : duohmm.cpp
// Author      : Jared O'Connell
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "hapmodule.h"
#include "pedigree.h"
#include "hmm.h"
#include "pedhap.h"
#include <boost/progress.hpp>

using namespace std;

int main(int argc,char **argv) {

  if(argc!=4) {
    cout << "Usage: duohmm hapfile famfile genetic_map --corrected output --recombination output --genotypingError output" << endl;
    return(0);
  }

  cout << "Reading haplotypes from " << argv[1] <<".haps" << endl;
  Haplotypes haps(argv[1]);

  cout << "Reading pedigree information from " << argv[2] << endl;

  pedigree p(argv[2],haps.ids);
  
  vector<string> idorder;
  p.orderSamples(p.pedigrees[0],idorder);
  for(int i=0;i<idorder.size();i++)
    cout << i <<": "<< idorder[i]<<endl;

  geneticMap gm(argv[3]);

  vector<double> cM;
  gm.interpolate(haps.positions,cM);
  for(int i=0;i<30;i++)
    cout << haps.positions[i] << "\t" << cM[i] << endl;

  DuoHMM duo(haps.positions,gm);
  cout << "RHO:"<<endl;
  for(int i=0;i<20;i++) {
    cout << duo.male_rho[i]<<"\t" << duo.male_norho[i]<<"\t" << duo.female_rho[i]<<"\t" << duo.female_norho[i]<<endl;
  }

  int i1=43;
  int i2=46;
  cout << haps.ids[i1] << "\t" << haps.ids[i2] << endl;
  //checks haplotypes are initialsed
  for(int i=0;i<30;i++)
    cout << (int)haps.H[i1*2][i] << " " << (int)haps.H[i1*2+1][i] << " " << (int)haps.H[i2*2][i] << " " << (int)haps.H[i2*2+1][i] << endl;
  cout << endl;

  int start = 4850;
  for(int i=start;i<start+20;i++)
    cout << (int)haps.H[i1*2][i] << " " << (int)haps.H[i1*2+1][i] << " " << (int)haps.H[i2*2][i] << " " << (int)haps.H[i2*2+1][i] << endl;
  cout << haps.nsnp << " " << haps.positions.size()<<endl;

  cout << endl;
  
  cout << "Setting haps..." <<endl;
  duo.setHaps(haps.getHap("ID178"),haps.getHap("ID276"),"2");  
  cout << "EM"<<endl;  
  duo.EM(10);
  cout << "viterbi" <<endl;
  duo.viterbi();

  /*
  for(int i=0;i<20;i++) {
    cout << duo.scale[i]<<endl;
    for(int j=0;j<4;j++) 
      cout << duo.alpha[j][i]<<"\t";
    cout << endl;
    for(int j=0;j<4;j++) 
      cout << duo.beta[j][i]<<"\t";
    cout << endl;
    for(int j=0;j<4;j++) 
      cout << duo.posterior[j][i]<<"\t";
    cout << endl;
  }
  cout << endl;
  for(int i=haps.nsnp-20;i<haps.nsnp;i++) {
    cout << duo.scale[i]<<endl;
    for(int j=0;j<4;j++) 
      cout << duo.alpha[j][i]<<"\t";
    cout << endl;
    for(int j=0;j<4;j++) 
      cout << duo.beta[j][i]<<"\t";
    cout << endl;
    for(int j=0;j<4;j++) 
      cout << duo.posterior[j][i]<<"\t";
    cout << endl;
  }
  */
  cout << duo.rho[100] << " " << duo.rho[1000] << endl;
  cout << endl;
  start = 4850;
  for(int i=start;i<start+20;i++) {
    cout << i << "\t";
    for(int j=0;j<4;j++) 
      cout << duo.posterior[j][i]<<"\t";
    cout << "\t" << duo.stateseq[i] << endl;
  }

  cout << "Making pedhap"<<endl;
  pedhap ph(argv[1],argv[2],argv[3]);
  ph.correct();

  duo.setHaps(ph.haps->getHap("ID178"),ph.haps->getHap("ID276"),"2");  
  cout << "EM"<<endl;  
  duo.EM(10);
  cout << "viterbi" <<endl;
  duo.viterbi();
  start = 4850;
  for(int i=start;i<start+20;i++) {
    cout << i << "\t";
    for(int j=0;j<4;j++) 
      cout << duo.posterior[j][i]<<"\t";
    cout << "\t" << duo.stateseq[i] << endl;
  }
  ph.haps->writeHaps(argv[1]);
  cout <<"\nAll tests completed.\nExiting..."<<endl;  return(0);
  
}

