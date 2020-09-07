//$Id$

#include "hmm.h"

geneticMap::geneticMap(string fname){
  pos.push_back(0);
  cM.push_back(0.0);    
  if(fname.compare("")==0) {
    pos.push_back(4000e6);
    cM.push_back(4000. * 1.19);
  } else {
    ifstream inf2(fname.c_str());
    if(!inf2) {
      cerr << "Problem reading genetic map: " << fname << "\nExiting..."<<endl;
      exit(1);
    }
    int tmp1;
    float tmp2,tmp3;
    inf2.ignore(1000,'\n');
    while(inf2) {
      inf2 >> tmp1;
      inf2 >> tmp2;
      inf2 >> tmp3;    
      //    if(DEBUG>0) 
      //      cout << tmp1 <<"\t" <<tmp2<<"\t"<<tmp3<<endl;
      pos.push_back(tmp1);
      cM.push_back(tmp3);
    }
    float     maxpos = pos[pos.size()-1] + 5e6;
    float last_cm = cM[cM.size()-1] + 5*1.19;
    pos.push_back(maxpos);
    cM.push_back(last_cm);
    //    cout << maxpos << " " << last_cm<<endl;
  }
  nsnp = pos.size();
}

#ifdef SHAPEIT
geneticMap::geneticMap(genhap_set & GH) {
  pos.push_back(0);
  cM.push_back(0.0);
  for(size_t l=0;l<GH.mapG->vec_pos.size();l++) {
    pos.push_back(GH.mapG->vec_pos[l]->bp);
    cM.push_back(GH.mapG->vec_pos[l]->cm);
  }
  nsnp = pos.size();
}
#endif

geneticMap::geneticMap() {
  int maxpos = 1000000000;
  pos.push_back(0);
  pos.push_back(maxpos);
  cM.push_back(0.0);
  cM.push_back(maxpos/1000000. * 1.19);
  nsnp = 2;
}

float geneticMap::interpolate(int position) {
  int i=0;
  while(i<nsnp && pos[i]<=position) {
    i++;
  }
  assert(i>0);
  if(i>=nsnp) {
    i--;
#ifndef SHAPEIT
    cerr << "WARNING: marker at position "<<position<< " was outside the range of provided genetic map"<<endl;
    cerr << "Is your map the appropriate chromosome and build?" << endl;
#endif
  }

  return(  cM[i-1] + (cM[i]-cM[i-1])*(float)(position-pos[i-1])/(float)(pos[i]-pos[i-1]) );
}

int geneticMap::interpolate(vector<int> & positions,vector<float> & output){
  output.resize(positions.size());

  for(size_t i=0;i<positions.size();i++) {
    output[i] = interpolate(positions[i]);
// #ifdef DEBUG
//     cerr<<positions[i]<<"\t"<<output[i]<<endl;
// #endif
  }

  return(0);
}

void DuoHMM::setIterations(int n) {
  NITERATION = n;
}

DuoHMM::DuoHMM(vector<int> & positions, geneticMap & gm)
{
  NITERATION=100;
  float r;
  male_multiplier = 0.7539868;
  female_multiplier = 1.2460132; 
  K=4;
  nsnp=positions.size();
  cM.resize(nsnp);
  gm.interpolate(positions,cM);
  male_length = (cM[nsnp-1]* male_multiplier)/100.;
  female_length = (cM[nsnp-1] * female_multiplier)/100.;

  male_rho.resize(nsnp);
  male_norho.resize(nsnp);
  female_rho.resize(nsnp);
  female_norho.resize(nsnp);

  for(int i=0;i<nsnp-1;i++) {

    r = (cM[i+1]-cM[i])/100.;
    if(r<=0.0) r = 1e-13;//hack to handle positions with same genetic location (which shouldnt be in the snps in the first place)
    male_norho[i] = exp(-r * male_multiplier);
    male_rho[i] = 1. - male_norho[i];
    female_norho[i] = exp(-r * female_multiplier);
    female_rho[i] = 1. - female_norho[i];
  }

  scale.resize(nsnp);

  alpha.assign(K,vector<float>(nsnp));
  beta.assign(K,vector<float>(nsnp));
  posterior.assign(K,vector<float>(nsnp));
  
  // alpha.resize(K);
  // beta.resize(K);
  // posterior.resize(K);
  // for(int i=0;i<K;i++) {
  //   alpha[i].resize(nsnp);
  //   beta[i].resize(nsnp);
  //   posterior[i].resize(nsnp);
  // }

  trans_posterior.resize(K);
  for(int i=0;i<K;i++)  {
    trans_posterior[i].resize(K);
    for(int j=0;j<K;j++) 
      trans_posterior[i][j].resize(nsnp,0.0);
  }

  error = 0.001;
  switch1 = 0.0005;
  switch2 = 0.0005;
}

  
int DuoHMM::forward() {
  int l = 0;
  float noerror = 1. - error;
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
  float tmp;
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
  float noerror = 1. - error;
  int l = nsnp-1;
  for(int i=0;i<K;i++) beta[i].assign(nsnp,0.0);

  //initial step
  for(int i=0;i<K;i++) beta[i][l] = scale[l];

  //induction
  
  float tmp;
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

int DuoHMM::setHaps(vector<bool> *parental_haplotypes,vector<bool> *child_haplotypes,string sex) {
  assert(parental_haplotypes!=NULL && child_haplotypes!=NULL);
  //  cout << sex << endl;
  assert(sex.compare("1")==0 || sex.compare("2")==0);
  parent = parental_haplotypes;
  child = child_haplotypes;
  if(sex.compare("1")==0)  {
    rho = &male_rho[0];
    genetic_length = male_length;
    multi = male_multiplier;
  }
  else {
    rho = &female_rho[0];
    genetic_length = female_length;
    multi = female_multiplier;
  }
  error = 0.001;
  switch1 = 0.0005;
  switch2 = 0.0005;
  
  return(0);
}

int DuoHMM::EM(int niteration) {


  float switch1_old=switch1;
  float switch2_old=switch2;
  float error_old=error;
  float tol=0.00001;
  float dif=2*tol;
  int i=0;
  while(i<niteration && dif>tol) {
    if(DEBUG>1) {
      cout << "ITERATION " <<i<<endl;
      cout << "error = "<< error << endl;
      cout << "switch1 = "<< switch1 << endl;
      cout << "switch2 = "<< switch2 << endl;
    }

    estep();
    mstep();

    dif = abs(switch1_old-switch1)+abs( switch2_old-switch2)+abs(error_old-error);
    switch1_old=switch1;
    switch2_old=switch2;
    error_old=error;
    i++;
  }

  return(0);
}

int DuoHMM::mstep() {

  float ngenerror = 0.0;
  float nswitch1 = 0.0;
  float nswitch2 = 0.0;
  int nhet1=0;
  int nhet2=0;
  int nhet3=0;
  genError.assign(nsnp,0.0);
  for(int l=0;l<nsnp;l++) {
    if(parent[0][l]!=parent[1][l]) nhet1++;
    if(child[0][l]!=child[1][l]) nhet2++;
    if(child[0][l]!=child[1][l] || parent[0][l]!=parent[1][l]) nhet3++;

    //genotyping error counts
    for(int i=0;i<K/2;i++) {
      for(int j=0;j<K/2;j++) {
	int  idx = i*2+j;
	if(child[j][l]!=parent[i][l])  {
	  ngenerror +=  posterior[idx][l];
	  genError[l] += posterior[idx][l];
	}
      }
    }
  }

  for(int l=0;l<nsnp-1;l++) {    //switch error counts
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

  //since these events are only detectable at hetorzygote sites, mstep denominator should probably be nhets rather than nsnp
  //in practice makes little difference

  /*
    switch1 =  (nswitch1 - genetic_length) / (nhet1 - 2*genetic_length);
    if(switch1<0.0) switch1 = 0.0;
    switch2 = nswitch2 / nhet2;
    error = ngenerror / nhet3;
  */

  switch1 =  (nswitch1 - genetic_length) / (nsnp - 2*genetic_length);
  if(switch1<0.0) switch1 = 0.0;
  switch2 = nswitch2 / nsnp;
  error = ngenerror / nsnp;

  return(0);
}

int DuoHMM::estep() {
  forward();
  backward();
  float noerror = 1. - error;

  float denominator;
  for(int l=0;l<nsnp;l++) {
    denominator = 0.0;
    for(int i=0;i<K;i++) {
      posterior[i][l] = alpha[i][l] * beta[i][l];
      denominator+=posterior[i][l];
    }
    for(int i=0;i<K;i++) posterior[i][l] /= denominator;
  }

  float aij;
  for(int l=0;l<(nsnp-1);l++) {
    assert(l<(nsnp-1));
    if(!(rho[l]>=0.0 && rho[l]<1.)) {
      cerr << "ERROR BAD RECOMBINATION RATE AT MARKER " << l << " rho = "<<rho[l]<<endl;
      cerr << cM[l] << " " << cM[l+1]<<endl;
      exit(1);
    }

    denominator = 0.0;

    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {
	int idx2 = i2*2+j2;

	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    int idx1 = i1*2+j1;
	    trans_posterior[idx1][idx2][l] = alpha[idx1][l]*beta[idx2][l+1];

	    if(i1==i2) 
	      aij =  (1-rho[l])*(1-switch1) + rho[l]*switch1 ; 
	    else 
	      aij = rho[l]*(1-switch1) + (1. - rho[l])*switch1 ; 

	    if(j1==j2) 
	      aij *= (1. - switch2);
	    else 
	      aij *= switch2;

	    trans_posterior[idx1][idx2][l] *= aij;

	    //emission
	    if(child[j2][l+1]==parent[i2][l+1])  
	      trans_posterior[idx1][idx2][l] *= noerror;
	    else 
	      trans_posterior[idx1][idx2][l] *= error;

	    denominator +=  trans_posterior[idx1][idx2][l];
	  }
	}

      }
    }

    for(int i2=0;i2<2;i2++) {
      for(int j2=0;j2<2;j2++) {
	int idx2 = i2*2+j2;
	for(int i1=0;i1<2;i1++) {//transition probs
	  for(int j1=0;j1<2;j1++) {
	    int idx1 = i1*2+j1;
	    //	    float tmp = trans_posterior[idx1][idx2][l];
	    trans_posterior[idx1][idx2][l] /= denominator;
	    //	    if(trans_posterior[idx1][idx2][l] > 0.5 && i1!=i2) 
	    //	      cout << l << "\t" << idx1 << "\t" << idx2 << "\t" << tmp << "\t" << denominator << endl;	    
	  }
	}
      }
    }

  }

  return(0);
}


int DuoHMM::viterbi() {
  int l = 0;
  float logerror = log(error);
  float lognoerror = log(1. - error);
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
  vector<float> tmp(4,0.0);
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

  if(DEBUG>1) {
    for(int l=0;l<nsnp;l++) 
      cout << l <<  "\t" << (int)stateseq[l] << endl;
  }
  return(0);
}

int DuoHMM::estimateRecombination() {

  recombinationMap.assign(nsnp,0.0);
  EM(NITERATION);
  float noerror = 1. - error;
  float r,rho2;
  int prevhet=-1;
  int l=0;
  while(parent[0][l]==parent[1][l]) l++;
  prevhet=l;

  vector <float> p;

  while(l<nsnp) {

    float recp=0;
    float denominator=0,p1,p2;
    while(parent[0][l]==parent[1][l] && l<nsnp) l++;

    if(l<nsnp) {
      r = multi * (cM[l]-cM[prevhet])/100.;    
      if(r<=0.0) r = 1e-13;//hack to handle positions with same genetic location (which shouldnt be in the snps in the first place)
      rho2 = 1 - exp(-r);

      for(int i2=0;i2<2;i2++) {
	for(int j2=0;j2<2;j2++) {

	  int idx2 = i2*2+j2;

	  for(int i1=0;i1<2;i1++) {//transition probs
	    for(int j1=0;j1<2;j1++) {

	      int idx1 = i1*2+j1;
	      if(child[j2][l]==parent[i2][l])  
		p1 = noerror*alpha[idx1][prevhet]*beta[idx2][l];
	      else
		p1 = error*alpha[idx1][prevhet]*beta[idx2][l];
	      p2 = p1;
	      
	      if(i1==i2) {
		p1 *= rho2*switch1; 
		p2 *= ( (1-rho2)*(1-switch1) + rho2*switch1 ); 
	      }
	      else {
		p1 *= rho2*(1-switch1);
		p2 *= ( rho2*(1-switch1) + (1. - rho2)*switch1 ); 
	      }
	      if(j1==j2)  {
		p1 *= (1. - switch2);
		p2 *= (1. - switch2);
	      }
	      else { 
		p1 *= switch2;
		p2 *= switch2;
	      }
	      recp += p1;
	      denominator +=  p2;
	    }
	  }

	}
      }
      recp/=denominator;
    }
    else {
      recp = 0.0;
    }
    if(!(recp>=0.0 && recp<=1.0)){
      cout << "ERROR:Invalid probability encountered"<<endl<<recp<<endl;
      cout << l << endl;
      cout <<  cM[l]<<" "<<cM[prevhet]<< endl;
      cout << r << endl;
      cout << rho2 << endl;
      cout << switch1 << endl;
      cout << switch2 << endl;
      exit(1);
    };

    for(int i=prevhet;i<l;i++) recombinationMap[i] = recp;
    prevhet=l;
    l++;
  }
  return 0;
}
