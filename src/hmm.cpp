#include "hmm.h"

#define DEBUG 0

geneticMap::geneticMap(string fname){
  io::filtering_istream inf2;
  ifstream blah(fname.c_str(),ios_base::in);
  inf2.push(blah);
  int tmp1;
  double tmp2,tmp3;
  inf2.ignore(1000,'\n');
  nsnp = 0;
  while(inf2) {
    inf2 >> tmp1;
    inf2 >> tmp2;
    inf2 >> tmp3;    
    //    if(DEBUG>0) 
    //    cout << tmp1 <<"\t" <<tmp2<<"\t"<<tmp3<<endl;
    pos.push_back(tmp1);
    cM.push_back(tmp3);
    nsnp++;
  }
}

double geneticMap::interpolate(int position) {
  int i=0;
  while(pos[i]<position) i++;
  return(  cM[i-1] + (cM[i]-cM[i-1])*(double)(position -pos[i-1])/(double)(pos[i]-pos[i-1]) );
}

int geneticMap::interpolate(vector<int> & positions,vector<double> & output){
  output.resize(positions.size());

  for(int i=0;i<positions.size();i++) 
    output[i] = interpolate(positions[i]);

  return(0);
}

DuoHMM::DuoHMM(vector<int> & positions, geneticMap & gm)
{
  double r;
  male_multiplier = 0.7539868;
  female_multiplier = 1.2460132; 
  K=4;
  nsnp=positions.size();
  cM.resize(nsnp);
  gm.interpolate(positions,cM);
  genetic_length = cM[nsnp-1];
  male_rho.resize(nsnp);
  male_norho.resize(nsnp);
  female_rho.resize(nsnp);
  female_norho.resize(nsnp);

  for(int i=0;i<nsnp-1;i++) {
    r = (cM[i+1]-cM[i])/100.;
    male_norho[i] = exp(-r * male_multiplier);
    male_rho[i] = 1. - male_norho[i];
    female_norho[i] = exp(-r * female_multiplier);
    female_rho[i] = 1. - female_norho[i];
  }

  scale.resize(nsnp);
  alpha.resize(K);
  beta.resize(K);
  posterior.resize(K);
  trans_posterior.resize(K);

  for(int i=0;i<K;i++)  {
    trans_posterior[i].resize(K);
    for(int j=0;j<K;j++) 
      trans_posterior[i][j].resize(nsnp,0.0);
  }

  for(int i=0;i<K;i++) {
    alpha[i].resize(nsnp);
    beta[i].resize(nsnp);
    posterior[i].resize(nsnp);
  }

  error = 0.001;
  switch1 = 0.0005;
  switch2 = 0.0005;
}

  
int DuoHMM::forward() {
  int l = 0;
  double noerror = 1. - error;
  scale.assign(nsnp,0.0);
  for(int i=0;i<K;i++) alpha[i].assign(nsnp,0.0);

  //initial step
  for(int i=0;i<K/2;i++) {
    for(int j=0;j<K/2;j++) {
      int  idx = i*2+j;
      if(child[j][l]==parent[i][l])  alpha[idx][l] = noerror;
      else alpha[idx][l] = error;
      alpha[idx][l] /= K;
      scale[l] += alpha[idx][l];
    }
  }
  scale[l] = 1./scale[l];
  for(int i=0;i<K;i++) alpha[i][l] *= scale[l];

  //induction
  int idx1,idx2;
  double tmp;
  for(l=1;l<nsnp;l++) {
    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {
	idx2 = i2*2+j2;
	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    idx1 = i1*2+j1;

	    if(i1==i2)  
	      tmp = (1-rho[l-1])*(1-switch1) + rho[l-1]*switch1; 
	    else 
	      tmp = rho[l-1]*(1-switch1) + (1. - rho[l-1])*switch1 ; 

	    if(j1==j2) 
	      tmp *= (1. - switch2);
	    else 
	      tmp *= switch2;			 

	    tmp *= alpha[idx1][l-1];
	    alpha[idx2][l]+=tmp;
	  }
	}

	//emission
	if(child[j2][l]==parent[i2][l])  
	  alpha[idx2][l] *= noerror;
	else 
	  alpha[idx2][l] *= error;

	scale[l] += alpha[idx2][l];	
      }
    }	    
    
    scale[l] = 1./scale[l];    
    for(int i=0;i<K;i++) alpha[i][l] *= scale[l];
  }

  return(0);
}

int DuoHMM::backward() {
  double noerror = 1. - error;
  int l = nsnp-1;
  for(int i=0;i<K;i++) beta[i].assign(nsnp,0.0);

  //initial step
  for(int i=0;i<K;i++) beta[i][l] = scale[l];

  //induction
  
  double tmp;
  for(l=nsnp-2;l>=0;l--) {
    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {

	int idx2 = i2*2+j2;
	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    int idx1 = i1*2+j1;
	    if(i1==i2)  tmp = (1-rho[l])*(1-switch1) + rho[l]*switch1; 
	    else tmp = rho[l]*(1. -switch1) + (1. - rho[l])*switch1 ; 
	    if(j1==j2) tmp *= (1. - switch2);
	    else tmp *= switch2;			 
	    tmp *= beta[idx1][l+1];
	    if(child[j1][l+1]==parent[i1][l+1])  
	      tmp *= noerror;
	    else 
	      tmp *= error;
	    
	    beta[idx2][l]+=tmp;
	  }
	}

	beta[idx2][l] *= scale[l];
      }
    }	    
  }

  return(0);
}

int DuoHMM::setHaps(unsigned char **parental_haplotypes,unsigned char **child_haplotypes,string sex) { 
  //  cout << sex << endl;
  assert(sex.compare("1")==0 || sex.compare("2")==0);
  parent = parental_haplotypes;
  child = child_haplotypes;
  if(sex.compare("1")==0) 
    rho = &male_rho[0];
  else 
    rho = &female_rho[0];
  error = 0.001;
  switch1 = 0.0005;
  switch2 = 0.0005;

  return(0);
}

int DuoHMM::EM(int niteration) {

  for(int i=0;i<niteration;i++) {
    if(DEBUG>1) {
      cout << "error = "<< error << endl;
      cout << "switch1 = "<< switch1 << endl;
      cout << "switch2 = "<< switch2 << endl;
    }
    estep();
    mstep();
  }

  return(0);
}

int DuoHMM::mstep() {

  double ngenerror = 0.0;
  double nswitch1 = 0.0;
  double nswitch2 = 0.0;
  int nhet1=0;
  int nhet2=0;
  int nhet3=0;

  for(int l=0;l<nsnp;l++) {
    if(parent[0][l]!=parent[1][l]) nhet1++;
    if(child[0][l]!=child[1][l]) nhet2++;
    if(child[0][l]!=child[1][l] || parent[0][l]!=parent[1][l]) nhet3++;

    //genotyping error counts
    for(int i=0;i<K/2;i++) {
      for(int j=0;j<K/2;j++) {
	int  idx = i*2+j;
	if(child[j][l]!=parent[i][l])  
	  ngenerror += posterior[idx][l];
      }
    }

    //switch error counts
    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {
	int idx2 = i2*2+j2;
	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    int idx1 = i1*2+j1;

	    if(i1!=i2) 
	      nswitch1 += trans_posterior[idx1][idx2][l];
	    if(j1!=j2)
	      nswitch2 += trans_posterior[idx1][idx2][l];

	  }
	}
      }
    }
  }

  if(DEBUG>1) {
    cout << "nswitch1 = " << nswitch1<<endl;
    cout << "nswitch2 = " << nswitch2<<endl;
    cout << "ngenerror = " << ngenerror<<endl;
  }

  switch1 =  (nswitch1 - genetic_length) / (nhet1 - 2*genetic_length);
  if(switch1<0.0) switch1 = 0.0;
  switch2 = nswitch2 / nhet2;

  switch1 = nswitch1 / nhet1;
  switch2 = nswitch2 / nhet2;
  error = ngenerror / nhet3;

  return(0);
}

int DuoHMM::estep() {
  forward();
  backward();
  double noerror = 1. - error;

  double denominator;
  for(int l=0;l<nsnp;l++) {
    denominator = 0.0;
    for(int i=0;i<K;i++) {
      posterior[i][l] = alpha[i][l] * beta[i][l];
      denominator+=posterior[i][l];
    }
    for(int i=0;i<K;i++) posterior[i][l] /= denominator;
  }

  double aij;
  for(int l=0;l<nsnp-1;l++) {
    denominator = 0.0;
    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {
	int idx2 = i2*2+j2;
	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    int idx1 = i1*2+j1;
	    trans_posterior[idx1][idx2][l] = alpha[idx1][l]*beta[idx2][l];
	    if(i1==i2) 
	      aij = ( (1-rho[l-1])*(1-switch1) + rho[l-1]*switch1 ); 
	    else 
	      aij = rho[l-1]*(1-switch1) + (1. - rho[l-1])*switch1 ; 
	    if(j1==j2) 
	      aij *= (1. - switch2);
	    else 
	      aij *= switch2;			 
	    trans_posterior[idx1][idx2][l] *= aij;
	    //emission
	    if(child[j2][l]==parent[i2][l])  
	      trans_posterior[idx1][idx2][l] *= noerror;
	    else 
	      trans_posterior[idx1][idx2][l] *= error;
	    denominator +=  trans_posterior[idx1][idx2][l];
	  }
	}
      }
    }
    for(int i=0;i<K;i++)
      for(int j=0;j<K;j++)
	trans_posterior[i][j][l] /= denominator;
  }

  return(0);
}


int DuoHMM::viterbi() {
  int l = 0;
  double logerror = log(error);
  double lognoerror = log(1. - error);
  stateseq.resize(nsnp);
  backtrack.resize(K);
  for(int i=0;i<K;i++) backtrack[i].assign(nsnp,0);

  //initial step
  for(int i=0;i<K/2;i++) {
    for(int j=0;j<K/2;j++) {
      int  idx = i*2+j;
      if(child[j][l]==parent[i][l])  alpha[idx][l] = lognoerror;
      else alpha[idx][l] = logerror;
      alpha[idx][l] -= log(K);
    }
  }

  //induction
  int idx1,idx2;
  vector<double> tmp(4,0.0);
  for(l=1;l<nsnp;l++) {
    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {
	idx2 = i2*2+j2;
	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    idx1 = i1*2+j1;

	    if(i1==i2)  
	      tmp[idx1] = log( (1-rho[l-1])*(1-switch1) + rho[l-1]*switch1 ); 
	    else 
	      tmp[idx1] = log(rho[l-1]*(1-switch1) + (1. - rho[l-1])*switch1); 

	    if(j1==j2) 
	      tmp[idx1] += log(1. - switch2);
	    else 
	      tmp[idx1] += log(switch2);

	    tmp[idx1] += alpha[idx1][l-1];
	  }
	}
	int maxind = argmax(tmp);
	backtrack[idx2][l]=maxind;
	//emission
	if(child[j2][l]==parent[i2][l])  
	  alpha[idx2][l] = tmp[maxind]+lognoerror;
	else 
	  alpha[idx2][l] = tmp[maxind]+logerror;

	if(DEBUG>1) {
	  cout.precision(4);
	  if(l>3787&&l<3797) {
	    cout << l << "\t";
	    for(int i=0;i<4;i++)
	      cout << std::fixed <<alpha[i][l] << "\t";
	    cout << endl;
	  }
	}

      }
    }	    
  }

  //backtrack
  int maxind=0;
  for(int i=1;i<K;i++) 
    if(alpha[i][nsnp-1]>alpha[maxind][nsnp-1])
      maxind=i;
  stateseq[nsnp-1]=maxind;

  for(l=nsnp-2;l>=0;l--) 
    stateseq[l] = backtrack[stateseq[l+1]][l+1];

  return(0);
}

