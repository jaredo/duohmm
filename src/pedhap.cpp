#include "pedhap.h"

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
  trio =  new TrioHMM(haps->positions,*gm);
  nsnp = haps->positions.size();
}

//constructor
pedhap::pedhap(string hap_filename,string pedigree_filename) {
  haps =  new Haplotypes(hap_filename);
  ped = new pedigree(pedigree_filename,haps->ids);
  gm = new geneticMap();
  duo =  new DuoHMM(haps->positions,*gm);
  trio =  new TrioHMM(haps->positions,*gm);
  nsnp = haps->positions.size();
}

int pedhap::phase(string child) {
  if(DEBUG>0)  cout << child << endl;
  assert(ped->sampleinfo.count(child));
  string dad = ped->sampleinfo[child].dad;
  string mum = ped->sampleinfo[child].mum;
  assert(ped->sampleinfo.count(dad)  ||  ped->sampleinfo.count(mum));
  unsigned char **d,**m,**p;
  unsigned char **c = haps->getHap(child);
  cout << "Phasing " << child << " - " << dad << " - " << mum << endl;
  pair<string,string> key,key1,key2;
  if(dad.compare("0")!=0 && mum.compare("0")!=0) {
    d = haps->getHap(dad);
    m = haps->getHap(mum);
    //    cout << (int)d[0][0] << (int)d[0][1]<< (int)m[0][0]<< (int)m[0][1]<< (int)c[0][0] <<  (int)c[0][1] << endl;
    trio->setHaps(d,m,c);
    trio->EM(NITERATION);
    trio->viterbi();
    key1 = make_pair(dad,child);
    key2 = make_pair(mum,child);
    stateseq[key1].resize(nsnp);
    stateseq[key2].resize(nsnp);  
    
    for(int l=0;l<nsnp;l++) {
      unsigned char s = trio->stateseq[l];
      stateseq[key1][l] = 2*(s/4) + s%2;
      stateseq[key2][l] = 2*((s%4)/2) + ((s%2)+1)%2;
    }
  } 
  else {
    if(dad.compare("0")!=0) {
      p = haps->getHap(dad);  
      duo->setHaps(p,c,"1");
      key = make_pair(dad,child);
    }
    if(mum.compare("0")!=0) {
      p = haps->getHap(mum);  
      duo->setHaps(p,c,"2");
      key = make_pair(mum,child);
    }
    duo->EM(NITERATION);
    duo->viterbi();    
    stateseq[key] = duo->stateseq;
  } 
  
  if(dad.compare("0")!=0) {
    key1 = make_pair(dad,child);
    int haporder = 0;
    for(int l=0;l<nsnp;l++) 
      if(stateseq[key1][l]%2!=haporder && c[0][l]!=c[1][l]) 
	unswitch(c,l,l+1);
  } else if(mum.compare("0")!=0) {
    key2 = make_pair(mum,child);
    int haporder = 0;
    for(int l=0;l<nsnp;l++) 
      if(stateseq[key2][l]%2!=haporder && c[0][l]!=c[1][l]) 
	unswitch(c,l,l+1);
  }

  if(dad.compare("0")!=0 && mum.compare("0")!=0) {
    trio->setHaps(d,m,c);
    trio->EM(NITERATION);
    trio->viterbi();
    for(int l=0;l<nsnp;l++) {
      unsigned char s = trio->stateseq[l];
      stateseq[key1][l] = 2*(s/4) + s%2;
      stateseq[key2][l] = 2*((s%4)/2) + ((s%2)+1)%2;
    }
  } 
  else {
    duo->EM(NITERATION);
    duo->viterbi();    
    stateseq[key] = duo->stateseq;
  } 

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

  /*
    for(set<string>::iterator it1=kids.begin();it1!=kids.end();it1++) {
    unsigned char **c = haps->getHap(*it1);
    duo->setHaps(p,c,ped->sampleinfo[parent].sex);  
    duo->EM(NITERATION);
    duo->viterbi();
    pair<string,string> key(parent,*it1);
    stateseq[key] = duo->stateseq;
    }
  */

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
      if(dad.compare("0")!=0 || mum.compare("0")!=0)  
	phase(*it2);
      /*
	if(dad.compare("0")!=0)  
	phase(dad,*it2);
	if(mum.compare("0")!=0)  
	phase(mum,*it2);
      */
    }

    for(vector<string>::iterator it2=idorder.begin();it2!=idorder.end();it2++) 
      minRecombinant(*it2);

  }

  return(0);
}


int pedhap::genotypingError(string outfile) {
  ofstream outf(outfile.c_str());
  cout << "Writing probable genotype errors to "<<outfile<<endl;
  outf << "ID FID MID RSID1 RSID2 POS PROB_ERROR\n";
  for(vector< set<string> >::iterator it1=ped->pedigrees.begin(); it1!=ped->pedigrees.end(); it1++) {

    vector<string> idorder;
    ped->orderSamples(*it1,idorder);

    for(vector<string>::iterator it2=idorder.begin();it2!=idorder.end();it2++) {
      string dad = ped->sampleinfo[*it2].dad;
      string mum = ped->sampleinfo[*it2].mum;
      cout << "Detecting errors for: " << *it2 << " - "<< dad << " - " << mum  << endl;

      if(dad.compare("0")!=0 && mum.compare("0")!=0) {

	unsigned char **c = haps->getHap(*it2);
	unsigned char **d = haps->getHap(dad);
	unsigned char **m = haps->getHap(mum);
	trio->setHaps(d,m,c);
	trio->EM(NITERATION);
	for(int l=0;l<nsnp;l++) {
	  double p = trio->genError[l];
	  if(p>0.1) {
	    outf << *it2 << " " << dad << " " << mum;
	    outf << " " << haps->rsid2[l]<< " " << haps->positions[l];
	    outf << " " << trio->genError[l] << endl;
	  }
	}

      } else if(dad.compare("0")!=0) {

	unsigned char **c = haps->getHap(*it2);
	unsigned char **d = haps->getHap(dad);
	duo->setHaps(d,c,"1");
	duo->EM(NITERATION);
	for(int l=0;l<nsnp;l++) {
	  double p = trio->genError[l];
	  if(p>0.1) {
	    outf << *it2 << " " << dad << " " << mum;
	    outf << " " << haps->rsid2[l]<< " " << haps->positions[l];
	    outf << " " << trio->genError[l] << endl;
	  }
	}
      } else if(mum.compare("0")!=0) {

	unsigned char **c = haps->getHap(*it2);
	unsigned char **m = haps->getHap(mum);
	duo->setHaps(m,c,"1");
	duo->EM(NITERATION);
	for(int l=0;l<nsnp;l++) {
	  double p = trio->genError[l];
	  if(p>0.1) {
	    outf << *it2 << " " << dad << " " << mum;
	    outf << " " << haps->rsid2[l]<< " " << haps->positions[l];
	    outf << " " << trio->genError[l] << endl;
	  }
	}
      }
    }
  }
  
  return 0;
}

int pedhap::recombinationMap(string outfile) { 
  ofstream outf(outfile.c_str());

  outf << "PARENT CHILD";
  for(int l=0;l<nsnp;l++) outf << " " << haps->positions[l];
  outf << endl;

  cout << "Calculating recombination map and writing to "<<outfile<<endl;
  for(vector< set<string> >::iterator it1=ped->pedigrees.begin(); it1!=ped->pedigrees.end(); it1++) {

    vector<string> idorder;
    ped->orderSamples(*it1,idorder);

    for(vector<string>::iterator it2=idorder.begin();it2!=idorder.end();it2++) {
	string dad = ped->sampleinfo[*it2].dad;
	string mum = ped->sampleinfo[*it2].mum;

      if(dad.compare("0")!=0 || mum.compare("0")!=0) {

	cout << "Mapping recombination for: " << *it2 << " - "<< dad << " - " << mum << endl;
	if(dad.compare("0")!=0 && mum.compare("0")!=0) {

	  unsigned char **c = haps->getHap(*it2);
	  unsigned char **d = haps->getHap(dad);
	  unsigned char **m = haps->getHap(mum);
	  trio->setHaps(d,m,c);
	  trio->estimateRecombination();
	  outf << *it2 << " " << dad << "\t";
	  for(int l=0;l<nsnp;l++) outf << trio->recombinationPat[l] << " ";
	  outf << endl << *it2 << " " << mum << "\t";
	  for(int l=0;l<nsnp;l++) outf << trio->recombinationMat[l] << " ";
	  outf << endl;

	} else if(dad.compare("0")!=0) {

	  unsigned char **c = haps->getHap(*it2);
	  unsigned char **d = haps->getHap(dad);
	  duo->setHaps(d,c,"1");
	  duo->estimateRecombination();
	  outf << *it2 << " " << dad<<"\t";
	  for(int l=0;l<nsnp;l++) outf << duo->recombinationMap[l] << " ";
	  outf << endl;
	} else if(mum.compare("0")!=0) {
	  unsigned char **c = haps->getHap(*it2);
	  unsigned char **m = haps->getHap(mum);
	  duo->setHaps(m,c,"2");
	  duo->estimateRecombination();
	  outf << *it2 << " " << mum<<"\t";
	  for(int l=0;l<nsnp;l++) outf << duo->recombinationMap[l] << " ";
	  outf << endl;
	}
      }
    }
  }  
  return 0;
}

