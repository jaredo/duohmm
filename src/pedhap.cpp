#include "pedhap.h"
#define DEBUG 5

int unswitch(unsigned char **h,int start,int stop) {
  for(int l=start;l<stop;l++)  {
    if(h[0][l]!=h[1][l]) {
      if(h[0][l]==1) {
	h[0][l]=0;
	h[1][l]=1;
      }
      else{
	h[0][l]=1;
	h[1][l]=0;
      }
    }
  }
  return(0);
}

//constructor
pedhap::pedhap(string hap_filename,string pedigree_filename,string gm_filename) {
  haps =  new Haplotypes(hap_filename);
  ped = new pedigree(pedigree_filename,haps->ids);
  gm = new geneticMap(gm_filename);
  duo =  new DuoHMM(haps->positions,*gm);
  nsnp = haps->positions.size();
}

int pedhap::phase(string parent,string child) {
  if(DEBUG>0)  cout << "Phasing "<<parent << " - "<<child << endl;
  assert(parent.compare("0")!=0);
  assert(ped->sampleinfo.count(parent) && ped->sampleinfo.count(child));

  unsigned char **p = haps->getHap(parent);
  unsigned char **c = haps->getHap(child);
  duo->setHaps(p,c,ped->sampleinfo[parent].sex);  
  duo->EM(10);
  duo->viterbi();
  
  int haporder = ped->sampleinfo[parent].sex.compare("1")!=0; //dad - 0     mum - 1
  for(int l=0;l<nsnp;l++) {
    if(duo->stateseq[l]%2!=haporder && c[0][l]!=c[1][l]) 
      unswitch(c,l,l+1);
  }

  if(DEBUG>2) {
    /*
    for(int l=0;l<nsnp;l++) {
      int maxind = 0;
      for(int i=1;i<4;i++) 
	if(duo->posterior[i][l]>duo->posterior[maxind][l]) 
	  maxind = i;
      
      if(duo->stateseq[l]!=maxind) {
	cout << l << "\t" << haps->positions[l] << "\t" << (int)duo->stateseq[l];
	for(int i=0;i<4;i++) 
	  cout << " " << duo->posterior[i][l];
	cout << endl;
      }
    }
    */
    cout << parent << " - " << child << endl;
    cout << "haporder = "<<haporder<<endl;
    cout << 0  << "\t"  << haps->positions[0] << "\t" << (int)duo->stateseq[0] << "\t" << (int)duo->stateseq[0]%2 << endl;
    for(int l=0;l<nsnp-1;l++) 
      if(duo->stateseq[l]!=duo->stateseq[l+1])  
	cout << l << "\t" << haps->positions[l] << "\t" << (int)duo->stateseq[l+1]  << "\t" << (int)duo->stateseq[l+1]%2  << endl;
    cout << nsnp-1   << "\t" << haps->positions[nsnp-1] << "\t" << (int)duo->stateseq[nsnp-1] << "\t" << (int)duo->stateseq[nsnp-1]%2 << endl;

  }

  duo->viterbi();
  //  pair<string,string> key(parent,child);
  //  stateseq[key] = duo->stateseq;

  return(0);
}

int pedhap::minRecombinant(string parent) {
  int minsib = 2;
  assert(ped->sampleinfo.count(parent));

  if(ped->sampleinfo[parent].kids.size()<minsib) 
    return(0);//not enough kids.

  if(ped->sampleinfo[parent].mum.compare("0")!=0 || ped->sampleinfo[parent].dad.compare("0")!=0)
    return(0);//has grand-parents = already phased

  if(DEBUG>0)  cout << "Finding minimum recombinant solution for "<<parent<<endl;
  unsigned char **p = haps->getHap(parent);
  set<string> kids = ped->sampleinfo[parent].kids;
  int prevhet=-1;


  for(set<string>::iterator it1=kids.begin();it1!=kids.end();it1++) {
    unsigned char **c = haps->getHap(*it1);
    duo->setHaps(p,c,ped->sampleinfo[parent].sex);  
    duo->EM(10);
    duo->viterbi();
    pair<string,string> key(parent,*it1);
    stateseq[key] = duo->stateseq;
  }

  for(int l=0;l<nsnp;l++) {
    if(p[0][l]!=p[1][l]){
      if(prevhet!=-1) {
	int nrec = 0;
	for(set<string>::iterator it1=kids.begin();it1!=kids.end();it1++) {
	  pair<string,string> key(parent,*it1);
	  if(stateseq[key][prevhet]/2 !=stateseq[key][l]/2)
	    nrec++;
	}
	if(nrec>minsib/2) {
	  if(DEBUG>0) cout << parent << " " << nrec <<  " " << l << " " << haps->positions[l]<<endl;
	  unswitch(p,l,nsnp);
	}
      }
      prevhet=l;
    }
  }

  return(0);
}

int pedhap::correct() {

  for(vector< set<string> >::iterator it1=ped->pedigrees.begin(); it1!=ped->pedigrees.end(); it1++) {

    vector<string> idorder;
    ped->orderSamples(*it1,idorder);

    for(vector<string>::iterator it2=idorder.begin();it2!=idorder.end();it2++) {
      string dad = ped->sampleinfo[*it2].dad;
      string mum = ped->sampleinfo[*it2].mum;
      if(dad.compare("0")!=0)  
	phase(dad,*it2);
      if(mum.compare("0")!=0)  
	phase(mum,*it2);
    }

    for(vector<string>::iterator it2=idorder.begin();it2!=idorder.end();it2++) 
      minRecombinant(*it2);

  }

  return(0);
}

