//$Id$

#include "hmm.h"

TrioHMM::TrioHMM(vector<int> & positions, geneticMap & gm)
{
    float r;
    NITERATION=100;
    K=8;
    male_multiplier = 0.7539868;
    female_multiplier = 1.2460132; 
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

    stateseq.resize(nsnp);
    scale.resize(nsnp);

    alpha.assign(K,vector<float>(nsnp));
    beta.assign(K,vector<float>(nsnp));
    posterior.assign(K,vector<float>(nsnp));
    // alpha.resize(K);
    // beta.resize(K);
    // posterior.resize(K);
    // for(int i=0;i<K;i++) {
    // 	alpha[i].resize(nsnp);
    // 	beta[i].resize(nsnp);
    // 	posterior[i].resize(nsnp);
    // }

    
    trans_posterior.resize(K);
    for(int i=0;i<K;i++)  {
	trans_posterior[i].resize(K);
	for(int j=0;j<K;j++) 
	    trans_posterior[i][j].resize(nsnp,0.0);
    }    
}

void TrioHMM::setIterations(int n) {
    NITERATION = n;
}

int TrioHMM::setHaps(unsigned char **dadptr,unsigned char **mumptr,unsigned char **childptr) {
    assert(dadptr!=nullptr && mumptr!=nullptr &&childptr!=nullptr);
    dad = dadptr;
    mum = mumptr;
    child = childptr;
    error = 0.005;
    switch_dad = 0.005;
    switch_mum = 0.005;
    switch_child = 0.005;

    return(0);
}

  
int TrioHMM::forward() {
    int l = 0;
    int kdad,kmum;
    float noerror = 1. - error;

    //initial step
    for(int i=0;i<2;i++) {
	for(int j=0;j<2;j++) {
	    for(int k=0;k<2;k++) {

		int  idx = i*4 + j*2 + k;
		if(k==0) {
		    kdad = 0;
		    kmum = 1;
		} else {
		    kdad = 1;
		    kmum = 0;	  
		}

		if(child[kdad][l]==dad[i][l] && child[kmum][l]==mum[j][l])
		    alpha[idx][l] = noerror;
		else 
		    alpha[idx][l] = error;

		alpha[idx][l] /= K;
	    }
	}
    }

    scale[0]=1.;

    //induction
    int idx1,idx2;
    float tmp;
    for(l=1;l<nsnp;l++) {
	scale[l]=0.0;

	for(int i2=0;i2<2;i2++) {//dad
	    for(int j2=0;j2<2;j2++) {//mum
		for(int k2=0;k2<2;k2++) {//kid
		    idx2 = i2*4 + j2*2 + k2;
		    alpha[idx2][l] = 0.0;
		    for(int i1=0;i1<2;i1++) {
			for(int j1=0;j1<2;j1++) {
			    for(int k1=0;k1<2;k1++) {
				idx1 = i1*4 + j1*2 + k1;

				if(i1==i2)  
				    tmp = (1-male_rho[l-1])*(1-switch_dad) + male_rho[l-1]*switch_dad; 
				else 
				    tmp = male_rho[l-1]*(1-switch_dad) + (1. - male_rho[l-1])*switch_dad ; 

				if(j1==j2)  
				    tmp *= (1-female_rho[l-1])*(1-switch_mum) + female_rho[l-1]*switch_mum; 
				else 
				    tmp *= female_rho[l-1]*(1-switch_mum) + (1. - female_rho[l-1])*switch_mum ; 

				if(k1==k2)
				    tmp *= (1. - switch_child);
				else 
				    tmp *= switch_child;

				tmp *= alpha[idx1][l-1];
				alpha[idx2][l]+=tmp;

			    }
			}
		    }
		    //emission
		    if(k2==0) {
			kdad = 0;
			kmum = 1;
		    } else {
			kdad = 1;
			kmum = 0;	  
		    }

		    if(child[kdad][l]==dad[i2][l] && child[kmum][l]==mum[j2][l])
			alpha[idx2][l] *= noerror;
		    else 
			alpha[idx2][l] *= error;

		    scale[l] += alpha[idx2][l];		

		}
	    }
	}
    
	scale[l] = 1./scale[l];    
	for(int i=0;i<K;i++) alpha[i][l] *= scale[l];
    }

    return(0);
}

int TrioHMM::backward() {
    int kdad,kmum;
    float noerror = 1. - error;

    for(int i=0;i<K;i++) beta[i].assign(nsnp,0.0);

    //initial step
    int l = nsnp-1;
    for(int i=0;i<K;i++) beta[i][l] = scale[l];

    //induction  
    float tmp;
    for(l=nsnp-2;l>=0;l--) {
	for(int i2=0;i2<2;i2++) {
	    for(int j2=0;j2<2;j2++) {
		for(int k2=0;k2<2;k2++) {

		    int idx2 = i2*4 + j2*2 + k2;
		    beta[idx2][l] = 0.0;
		    for(int i1=0;i1<2;i1++) {//transition probs
			for(int j1=0;j1<2;j1++) {
			    for(int k1=0;k1<2;k1++) {

				int idx1 = i1*4+j1*2+k1;
				if(i1==i2)  
				    tmp = (1-male_rho[l])*(1-switch_dad) + male_rho[l]*switch_dad; 
				else 
				    tmp = male_rho[l]*(1-switch_dad) + (1. - male_rho[l])*switch_dad ; 

				if(j1==j2)  
				    tmp *= (1-female_rho[l])*(1-switch_mum) + female_rho[l]*switch_mum; 
				else 
				    tmp *= female_rho[l]*(1-switch_mum) + (1. - female_rho[l])*switch_mum ; 

				if(k1==k2)
				    tmp *= (1. - switch_child);
				else 
				    tmp *= switch_child;
			 
				tmp *= beta[idx1][l+1];
				//emission
				if(k1==0) {
				    kdad = 0;
				    kmum = 1;
				} else {
				    kdad = 1;
				    kmum = 0;	  
				}
				if(child[kdad][l+1]==dad[i1][l+1] && child[kmum][l+1]==mum[j1][l+1])
				    tmp *= noerror;
				else 
				    tmp *= error;
		
				beta[idx2][l]+=tmp;
			    }
			}
		    }

		    beta[idx2][l]*=scale[l];

		}
	    }	    
	}
    }

    return(0);
}


int TrioHMM::EM(int niteration) {

    float switch_dad_old=switch_dad;
    float switch_mum_old=switch_mum;
    float switch_child_old=switch_child;
    float error_old=error;
    float tol=0.00001;
    float dif=2*tol;
    int i=0;

    while(i<niteration && dif>tol) {
	if(DEBUG>1) {
	    cout.precision(10);
	    cout << "\nITERATION "<<i<<endl;
	    cout << "error = "<< error << endl;
	    cout << "switch_dad = "<< switch_dad << endl;
	    cout << "switch_mum = "<< switch_mum << endl;
	    cout << "switch_child = "<< switch_child << endl;
	}

	estep();
	mstep();

	dif = abs(switch_dad_old-switch_dad)+abs( switch_mum_old-switch_mum)+abs(switch_child_old-switch_child)+abs(error_old-error);
	switch_dad_old=switch_dad;
	switch_mum_old=switch_mum;
	switch_child_old=switch_child;
	error_old=error;
	i++;
    }

    return(0);
}

int TrioHMM::mstep() {

    float ngenerror = 0.0;
    float nswitch_dad = 0.0;
    float nswitch_mum = 0.0;
    float nswitch_child = 0.0;
    int nhet_mum=0;
    int nhet_dad=0;
    int nhet_child=0;
    int nhet_inf=0;
    int kdad,kmum;
    float  nconsistent = 0;
    genError.assign(nsnp,0.0);

    for(int l=0;l<nsnp;l++) {

	if(dad[0][l]!=dad[1][l]) nhet_dad++;
	if(mum[0][l]!=mum[1][l]) nhet_mum++;
	if(child[0][l]!=child[1][l]) nhet_child++;
	if(child[0][l]!=child[1][l] || dad[0][l]!=dad[1][l] || mum[0][l]!=mum[1][l]) nhet_inf++;

	//genotyping error counts
	for(int i=0;i<2;i++) {
	    for(int j=0;j<2;j++) {
		for(int k=0;k<2;k++) {
		    int  idx = i*4 + j*2 + k;
		    //emission
		    if(k==0) {
			kdad = 0;
			kmum = 1;
		    } else {
			kdad = 1;
			kmum = 0;	  
		    }	  
		    if(child[kdad][l]!=dad[i][l] || child[kmum][l]!=mum[j][l] )  {
			ngenerror += posterior[idx][l];
			genError[l] += posterior[idx][l];
		    }
		    else
			nconsistent += posterior[idx][l];
		}
	    }
	}
    }

    for(int l=0;l<nsnp-1;l++) {
	//switch error counts
	for(int i2=0;i2<2;i2++) {
	    for(int j2=0;j2<2;j2++) {
		for(int k2=0;k2<2;k2++) {
		    int idx2 = i2*4 + j2*2 + k2;
		    for(int i1=0;i1<2;i1++) {//transition probs
			for(int j1=0;j1<2;j1++) {
			    for(int k1=0;k1<2;k1++) {
				int idx1 = i1*4 + j1*2 +k1;

				if(i1!=i2) 
				    nswitch_dad += trans_posterior[idx1][idx2][l];
				if(j1!=j2)
				    nswitch_mum += trans_posterior[idx1][idx2][l];
				if(k1!=k2)
				    nswitch_child += trans_posterior[idx1][idx2][l];
		
			    }
			}
		    }
		}
	    }
	}
    }

    if(DEBUG>1) {
	cout << "nswitch_mum = " << nswitch_mum    << " / " <<  nhet_mum <<endl;
	cout << "nswitch_dad = " << nswitch_dad     << " / " <<  nhet_dad <<endl;
	cout << "nswitch_child = " << nswitch_child    << " / " <<  nhet_child <<endl;
	cout << "ngenerror = " << ngenerror    << " / " <<  nhet_inf<<endl;
	cout << "nconsistent = " << nconsistent << endl;
	cout << male_length << " " << female_length << endl << endl;
    }

    nhet_dad=nsnp;
    nhet_child=nsnp;
    nhet_mum=nsnp;
    //since these events are only detectable at hetorzygote sites, mstep denominator should probably be nhets rather than nsnp
    //in practice makes little difference
    /*
      switch_dad =  (nswitch_dad - male_length) / (nhet_dad - 2*male_length);
      switch_mum =  (nswitch_mum - female_length) / (nhet_mum - 2*female_length);
      if(switch_dad<0.0) switch_dad = 0.0;
      if(switch_mum<0.0) switch_mum = 0.0;
      switch_child = nswitch_child / nhet_child;
      error = ngenerror / nhet_inf;
    */

    switch_dad =  (nswitch_dad - male_length) / (nsnp - 2*male_length);
    switch_mum =  (nswitch_mum - female_length) / (nsnp - 2*female_length);
    if(switch_dad<0.0) switch_dad = 0.0;
    if(switch_mum<0.0) switch_mum = 0.0;
    switch_child = nswitch_child / nsnp;
    error = ngenerror / nsnp;


    if(error<=0.0) error = 0.00001;

    return(0);
}

int TrioHMM::estep() {
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

    int kdad,kmum;
    float aij;
    for(int l=0;l<nsnp-1;l++) {
	denominator = 0.0;
	for(int i2=0;i2<2;i2++) {
	    for(int j2=0;j2<2;j2++) {
		for(int k2=0;k2<2;k2++) {
		    int idx2 = i2*4 + j2*2 +k2;
		    for(int i1=0;i1<2;i1++) {//transition probs
			for(int j1=0;j1<2;j1++) {
			    for(int k1=0;k1<2;k1++) {

				int idx1 = i1*4+j1*2+k1;
				trans_posterior[idx1][idx2][l] = alpha[idx1][l]*beta[idx2][l+1];
				if(i1==i2) 
				    aij = ( (1-male_rho[l])*(1-switch_dad) + male_rho[l]*switch_dad ); 
				else 
				    aij = male_rho[l]*(1-switch_dad) + (1. - male_rho[l])*switch_dad ; 

				if(j1==j2) 
				    aij *= ( (1-female_rho[l])*(1-switch_mum) + female_rho[l]*switch_mum ); 
				else 
				    aij *= female_rho[l]*(1-switch_mum) + (1. - female_rho[l])*switch_mum ; 

				if(k1==k2) 
				    aij *= (1. - switch_child);
				else 
				    aij *= switch_child;	
		 
				trans_posterior[idx1][idx2][l] *= aij;

				//emission
				if(k1==0) {
				    kdad = 0;
				    kmum = 1;
				} else {
				    kdad = 1;
				    kmum = 0;	  
				}

				if(child[kdad][l+1]==dad[i1][l+1] && child[kmum][l+1]==mum[j1][l+1])
				    trans_posterior[idx1][idx2][l] *= noerror;
				else 
				    trans_posterior[idx1][idx2][l] *= error;

				denominator +=  trans_posterior[idx1][idx2][l];

			    }
			}
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


int TrioHMM::viterbi() {
    int l = 0;
    float logerror = log(error);
    float lognoerror = log(1. - error);
    stateseq.resize(nsnp);
    backtrack.resize(K);
    for(int i=0;i<K;i++) backtrack[i].assign(nsnp,0);
    int kmum,kdad;

    //initial step
    for(int i=0;i<2;i++) {
	for(int j=0;j<2;j++) {
	    for(int k=0;k<2;k++) {
		int  idx = i*4 + j*2 + k;
		if(k==0) {
		    kdad = 0;
		    kmum = 1;
		} else {
		    kdad = 1;
		    kmum = 0;	  
		}
		if(child[kdad][l]==dad[i][l] && child[kmum][l]==mum[j][l])
		    alpha[idx][l] = lognoerror;
		else 
		    alpha[idx][l] = logerror;

		alpha[idx][l] -= K;
	    }
	}
    }


    //induction
    int idx1,idx2;
    vector<float> tmp(K,0.0);

    for(l=1;l<nsnp;l++) {

	for(int i2=0;i2<2;i2++) {//dad
	    for(int j2=0;j2<2;j2++) {//mum
		for(int k2=0;k2<2;k2++) {//kid

		    idx2 = i2*4 + j2*2 + k2;

		    for(int i1=0;i1<2;i1++) {
			for(int j1=0;j1<2;j1++) {
			    for(int k1=0;k1<2;k1++) {
				idx1 = i1*4 + j1*2 + k1;
				if(i1==i2)  
				    tmp[idx1] = log( (1-male_rho[l-1])*(1-switch_dad) + male_rho[l-1]*switch_dad ); 
				else 
				    tmp[idx1] = log( male_rho[l-1]*(1-switch_dad) + (1. - male_rho[l-1])*switch_dad ); 

				if(j1==j2)  
				    tmp[idx1] += log( (1-female_rho[l-1])*(1-switch_mum) + female_rho[l-1]*switch_mum ); 
				else 
				    tmp[idx1] += log( female_rho[l-1]*(1-switch_mum) + (1. - female_rho[l-1])*switch_mum ); 

				if(k1==k2)
				    tmp[idx1] += log(1. - switch_child);
				else 
				    tmp[idx1] += log(switch_child);
		      
				tmp[idx1] += alpha[idx1][l-1];
		
			    }	    
			}
		    }

		    int maxind = argmax(tmp);
		    backtrack[idx2][l]=maxind;
		    //emission
		    if(k2==0) {
			kdad = 0;
			kmum = 1;
		    } else {
			kdad = 1;
			kmum = 0;	  
		    }
		    if(child[kdad][l]==dad[i2][l] && child[kmum][l]==mum[j2][l])
			alpha[idx2][l] = tmp[maxind]+lognoerror;
		    else 
			alpha[idx2][l] = tmp[maxind]+logerror;
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


int TrioHMM::estimateRecombination()
{
    recombinationPat.assign(nsnp,0.0);
    recombinationMat.assign(nsnp,0.0);
    EM(NITERATION);
    estimateRecombinationDad();
    estimateRecombinationMum();
    return 0;
}

int TrioHMM::estimateRecombinationDad() {
    float r,rho_dad;
    float noerror = 1. - error;
    int prevhet=-1;
    int l=0;
    float p1,p2,recp;
    int kdad,kmum;
    while(dad[0][l]==dad[1][l]) l++;
    prevhet=l;

    while(l<nsnp) {
	l++;
	while(dad[0][l]==dad[1][l] && l<nsnp) l++;
	if(l<nsnp) {

	    r =  (cM[l]-cM[prevhet])/100.;    
	    if(r<=0.0) r = 1e-13;//hack to handle positions with same genetic location (which shouldnt be in the snps in the first place)

	    rho_dad = 1 - exp(-male_multiplier*r);
	    //      rho_dad = 0.5;
	    recp = 0.0;
	    float denominator=0;//recombination probability

	    for(int i2=0;i2<2;i2++) {
		for(int j2=0;j2<2;j2++) {
		    for(int k2=0;k2<2;k2++) {
			int idx2 = i2*4+j2*2+k2;

			for(int i1=0;i1<2;i1++) {//transition probs
			    for(int j1=0;j1<2;j1++) {
				for(int k1=0;k1<1;k1++) {
				    int idx1 = i1*4+j1*2+k1;

				    if(k2==0) {
					kdad = 0;
					kmum = 1;
				    } else {
					kdad = 1;
					kmum = 0;	  
				    }
				    p1 = alpha[idx1][prevhet]*beta[idx2][l];
				    if(child[kdad][l]==dad[i2][l] && child[kmum][l]==mum[j2][l])
					p1*=noerror;
				    else 
					p1*=error;

				    p2 = p1;

				    if(i1==i2) {
					p1 *= rho_dad*switch_dad ;
					p2 *= ( (1-rho_dad)*(1-switch_dad) + rho_dad*switch_dad ); 
				    }
				    else {
					p1 *= rho_dad*(1-switch_dad);
					p2 *= rho_dad*(1-switch_dad) + (1. - rho_dad)*switch_dad ; 
				    }
				    recp += p1;
				    denominator +=  p2;
				}
			    }
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
	    cout << rho_dad << endl;
	    cout << switch_dad << endl;
	    cout << switch_mum << endl;

	    exit(1);
	};
	for(int i=prevhet;i<l;i++) recombinationPat[i] = recp;
	prevhet = l;
    }
    return 0;
}

int TrioHMM::estimateRecombinationMum() {
    float r,rho_mum;
    float noerror = 1. - error;
    int prevhet=-1;
    int l=0;
    float p1,p2,recp;
    int kdad,kmum;

    while(mum[0][l]==mum[1][l]) l++;
    prevhet=l;
  
    while(l<nsnp) {
	l++;
	while(mum[0][l]==mum[1][l] && l<nsnp) l++;
	if(l<nsnp) {

	    r =  (cM[l]-cM[prevhet])/100.;    
	    if(r<=0.0) r = 1e-13;//hack to handle positions with same genetic location (which shouldnt be in the snps in the first place)
	    rho_mum = 1 - exp(-female_multiplier*r);
	    recp = 0.0;
	    float denominator=0;//recombination probability

	    for(int i2=0;i2<2;i2++) {
		for(int j2=0;j2<2;j2++) {
		    for(int k2=0;k2<2;k2++) {
			int idx2 = i2*4+j2*2+k2;

			for(int i1=0;i1<2;i1++) {//transition probs
			    for(int j1=0;j1<2;j1++) {
				for(int k1=0;k1<1;k1++) {
				    int idx1 = i1*4+j1*2+k1;

				    if(k2==0) {
					kdad = 0;
					kmum = 1;
				    } else {
					kdad = 1;
					kmum = 0;	  
				    }

				    if(child[kdad][l]==dad[i2][l] && child[kmum][l]==mum[j2][l])
					p1 = alpha[idx1][prevhet]*beta[idx2][l]*noerror;
				    else 
					p1 = alpha[idx1][prevhet]*beta[idx2][l]*error;

				    p2 = p1;

				    if(j1==j2) {
					p1 *= rho_mum*switch_mum;
					p2 *= ( (1-rho_mum)*(1-switch_mum) + rho_mum*switch_mum ); 
				    }
				    else {
					p1 *= rho_mum*(1-switch_mum);
					p2 *= rho_mum*(1-switch_mum) + (1. - rho_mum)*switch_mum ; 
				    }
				    recp += p1;
				    denominator +=  p2;
				}
			    }
			}

		    }
		}
	    }
	    recp/=denominator;    
	}
	else {
	    recp = 0.0;
	}
	assert(recp>=0.0 && recp<=1.0);
	for(int i=prevhet;i<l;i++) recombinationMat[i] = recp;
	prevhet = l;
    }

    return 0;
}
