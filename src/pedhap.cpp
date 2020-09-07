//$Id$

#include "pedhap.h"

int unswitch(vector<bool> *h, int start, int stop)
{
  for (int l = start; l < stop; l++)
  {
    if (h[0][l] != h[1][l])
    {
      if (h[0][l] == 1)
      {
        h[0][l] = 0;
        h[1][l] = 1;
      }
      else
      {
        h[0][l] = 1;
        h[1][l] = 0;
      }
    }
  }
  return (0);
}

#ifdef SHAPEIT
pedhap::pedhap(filter_writer &F, genhap_set &GH, string header1, string header2, int niteration)
{
  NITERATION = niteration;
  if (VERBOSE)
    cout << "converting haps" << endl;
  haps = new Haplotypes(F, GH);
  ped = new pedigree(F, GH, header1, header2, haps->ids);
  if (VERBOSE)
    cout << "reading genetic map." << endl;
  gm = new geneticMap(GH);
  duo = new DuoHMM(haps->positions, *gm);
  trio = new TrioHMM(haps->positions, *gm);
  nsnp = haps->positions.size();
}

#endif

//constructor
pedhap::pedhap(string hap_filename, string pedigree_filename, string gm_filename)
{
  haps = new Haplotypes(hap_filename);
  ped = new pedigree(pedigree_filename, haps->ids, 'f');
  gm = new geneticMap(gm_filename);
  duo = new DuoHMM(haps->positions, *gm);
  trio = new TrioHMM(haps->positions, *gm);
  nsnp = haps->positions.size();
}

//constructor
pedhap::pedhap(string hap_filename, string gm_filename, int niteration)
{
  haps = new Haplotypes(hap_filename);
  ped = new pedigree(hap_filename + ".sample", haps->ids, 's');
  gm = new geneticMap(gm_filename);
  duo = new DuoHMM(haps->positions, *gm);
  duo->setIterations(niteration);
  trio = new TrioHMM(haps->positions, *gm);
  trio->setIterations(niteration);
  nsnp = haps->positions.size();
  NITERATION = niteration;
}

pedhap::~pedhap()
{
  delete haps;
  delete ped;
  delete gm;
  delete duo;
  delete trio;
}

int pedhap::phase(string child)
{
  bool impute_missing_genotypes=false; //re-impute missing genotypes with a simple haploid imputation model
  if (DEBUG > 0)
    cout << child << endl;
  assert(ped->sampleinfo.count(child));
  string dad = ped->sampleinfo[child].dad;
  string mum = ped->sampleinfo[child].mum;
  string parent="";
  if (VERBOSE > 0)
  {
    cout << "Phasing " << child << " - " << dad << " - " << mum << endl;
  }
  assert(ped->sampleinfo.count(dad) || ped->sampleinfo.count(mum));
#ifdef SHAPEIT
  vector<bool> *kid_missing = haps->getMissing(child);
  vector<bool> *dad_missing = haps->getMissing(dad);
  vector<bool> *mum_missing = haps->getMissing(mum);
  vector<bool> *parent_missing = NULL;
  if (mum.compare("0") != 0)
    parent_missing = mum_missing;
  if (dad.compare("0") != 0)
    parent_missing = dad_missing;
#endif
  vector<bool> *d = NULL, *m = NULL, *p = NULL;
  vector<bool> *c = haps->getHap(child);
  pair<string, string> key, key1, key2;

  //trio
  if (dad.compare("0") != 0 && mum.compare("0") != 0)
  {
    d = haps->getHap(dad);
    m = haps->getHap(mum);
    //    cout << (int)d[0][0] << (int)d[0][1]<< (int)m[0][0]<< (int)m[0][1]<< (int)c[0][0] <<  (int)c[0][1] << endl;
    trio->setHaps(d, m, c);
    trio->EM(NITERATION);
    trio->viterbi();
    key1 = make_pair(dad, child);
    key2 = make_pair(mum, child);
    stateseq[key1].resize(nsnp);
    stateseq[key2].resize(nsnp);

    for (int l = 0; l < nsnp; l++)
    {
      unsigned char s = trio->stateseq[l];
      stateseq[key1][l] = 2 * (s / 4) + s % 2;
      stateseq[key2][l] = 2 * ((s % 4) / 2) + ((s % 2) + 1) % 2;
    }
  }
  else
  { //duos
    if (dad.compare("0") != 0)
    {
	parent=dad;
      p = haps->getHap(dad);
      duo->setHaps(p, c, "1");
      key = make_pair(dad, child);
    }
    if (mum.compare("0") != 0)
    {
	parent=mum;
      p = haps->getHap(mum);
      duo->setHaps(p, c, "2");
      key = make_pair(mum, child);
    }
    duo->EM(NITERATION);
    duo->viterbi();
    stateseq[key] = duo->stateseq;
  }

  if (dad.compare("0") != 0)
  {
    key1 = make_pair(dad, child);
    int haporder = 0;
    for (int l = 0; l < nsnp; l++)
      if (stateseq[key1][l] % 2 != haporder && c[0][l] != c[1][l])
        unswitch(c, l, l + 1);
  }
  else if (mum.compare("0") != 0)
  {
    key2 = make_pair(mum, child);
    int haporder = 1;
    for (int l = 0; l < nsnp; l++)
      if (stateseq[key2][l] % 2 != haporder && c[0][l] != c[1][l])
        unswitch(c, l, l + 1);
  }

  if (dad.compare("0") != 0 && mum.compare("0") != 0)
  {
    trio->setHaps(d, m, c);
    trio->EM(NITERATION);
    trio->viterbi();
    for (int l = 0; l < nsnp; l++)
    {
      unsigned char s = trio->stateseq[l];
      stateseq[key1][l] = 2 * (s / 4) + s % 2;
      stateseq[key2][l] = 2 * ((s % 4) / 2) + ((s % 2) + 1) % 2;

#ifdef SHAPEIT
      //re-impute missing genotypes
      //when entire trio/duo is missing. re-impute parents from population

      if(impute_missing_genotypes)
      {
	  if ((*kid_missing)[l] && (*dad_missing)[l])
	      haps->impute(dad,l);
	  if ((*kid_missing)[l] && (*mum_missing)[l])
	      haps->impute(mum,l);	  
      }
      
      bitset<3> b(s);
      int mum_src = b[1];
      int dad_src = b[2];	  
      //always impute a missing child directly from parents
      if ((*kid_missing)[l])
      {
        c[b[0]][l] = d[dad_src][l];
        c[(b[0] + 1) % 2][l] = m[mum_src][l];
        //	      assert(is_mendel_consistent(c[0][l]+c[1][l],d[0][l]+d[1][l],m[0][l]+m[1][l]));//debug
        //	      debugging
        //	      std::cerr << "MISSING KID "<< child <<" "<<dad<<" "<<mum<<" "<<haps->positions[l]<<" "<<(int)trio->stateseq[l]<<" "<<b[0]<<b[1]<<b[2];
        //	      std::cerr<<" "<<(int)c[0][l]<<" "<<(int)c[1][l]<<" "<<(int)m[0][l]<<" "<<(int)m[1][l]<<" "<<(int)d[0][l]<<" "<<(int)d[1][l]<<endl;
      }

      //if there is mendel inconsistency. impute missing parents from child
      if (!is_mendel_consistent(c[0][l] + c[1][l], d[0][l] + d[1][l], m[0][l] + m[1][l]))
      {
        if ((*dad_missing)[l] && (*mum_missing)[l])
        {
          d[dad_src][l] = c[b[0]][l];
          m[mum_src][l] = c[(b[0] + 1) % 2][l];
        }
        else if ((*dad_missing)[l])
        {
          d[dad_src][l] = c[b[0]][l];
          if (!is_mendel_consistent(c[0][l] + c[1][l], d[0][l] + d[1][l], m[0][l] + m[1][l]))
          {
            d[dad_src][l] = (d[dad_src][l] + 1) % 2;
          }
        }
        else if ((*mum_missing)[l])
        {
          m[mum_src][l] = c[(b[0] + 1) % 2][l];
          if (!is_mendel_consistent(c[0][l] + c[1][l], d[0][l] + d[1][l], m[0][l] + m[1][l]))
          {
            m[mum_src][l] = (m[mum_src][l] + 1) % 2;
          }
        }
      }
#endif
    }
  }
  else
  {
    duo->EM(NITERATION);
    duo->viterbi();
    stateseq[key] = duo->stateseq;
#ifdef SHAPEIT //re-impute missing genotypes
    for (int l = 0; l < nsnp; l++)
    {
	if ((*kid_missing)[l] && (*parent_missing)[l])
	{
	    haps->impute(parent,l);
	}
	unsigned char s = duo->stateseq[l];
	bitset<2> b(s);
	int src = b[1];
	int dst = b[0];
	if ((*kid_missing)[l])
	{
	    c[dst][l] = p[src][l];
	}
	if((*parent_missing)[l])
	{
	    p[src][l] = c[dst][l];
	}
    }
#endif
  }

  return (0);
}

int pedhap::minRecombinant(string parent)
{
  int minsib = 2;
  assert(ped->sampleinfo.count(parent));
  int nkid = ped->sampleinfo[parent].kids.size();
  if (nkid < minsib)
    return (0); //not enough kids.

  /*
    if(ped->sampleinfo[parent].mum.compare("0")!=0 || ped->sampleinfo[parent].dad.compare("0")!=0)
      return(0);//has grand-parents = already phased
    */

  if (DEBUG > 0)
    cout << "Finding minimum recombinant solution for " << parent << endl;
  vector<bool> *p = haps->getHap(parent);
  set<string> kids = ped->sampleinfo[parent].kids;
  int prevhet = -1;

  for (int l = 0; l < nsnp; l++)
  {
    if (p[0][l] != p[1][l])
    {
      if (prevhet != -1)
      {
        int nrec = 0;
        for (set<string>::iterator it1 = kids.begin(); it1 != kids.end(); it1++)
        {
          pair<string, string> key(parent, *it1);
          if (stateseq[key][prevhet] / 2 != stateseq[key][l] / 2)
            nrec++;
        }
        if (nrec > nkid / 2)
        {
          if (DEBUG > 0)
            cout << parent << " " << nrec << " " << l << " " << haps->positions[l] << endl;
          unswitch(p, l, nsnp);
        }
      }
      prevhet = l;
    }
  }

  return (0);
}

int pedhap::correct()
{
  int num_phaseable = 0, num_phased = 0;

//this just counts how many meioses we are going to phasew
#ifdef SHAPEIT
  for (vector<set<string>>::iterator it1 = ped->pedigrees.begin(); it1 != ped->pedigrees.end(); it1++)
  {
    vector<string> idorder;
    ped->orderSamples(*it1, idorder);
    for (vector<string>::iterator it2 = idorder.begin(); it2 != idorder.end(); it2++)
    {
      string dad = ped->sampleinfo[*it2].dad;
      string mum = ped->sampleinfo[*it2].mum;
      if (dad.compare("0") != 0 || mum.compare("0") != 0)
      {
        num_phaseable++;
      }
    }
  }
  cout << "\nCorrecting haplotypes based on pedigree structure [" << num_phased << "/" << num_phaseable << "]\r";
#endif

  //this does the phasing
  for (vector<set<string>>::iterator it1 = ped->pedigrees.begin(); it1 != ped->pedigrees.end(); it1++)
  {
    vector<string> idorder;
    ped->orderSamples(*it1, idorder);
    for (vector<string>::iterator it2 = idorder.begin(); it2 != idorder.end(); it2++)
    {
      string dad = ped->sampleinfo[*it2].dad;
      string mum = ped->sampleinfo[*it2].mum;
      if (dad.compare("0") != 0 || mum.compare("0") != 0)
      {
        num_phased++;
        phase(*it2);
#ifdef SHAPEIT
        cout << "Correcting haplotypes based on pedigree structure [" << num_phased << "/" << num_phaseable << "]\r";
#endif
      }
    }
    for (vector<string>::iterator it2 = idorder.begin(); it2 != idorder.end(); it2++)
    {
      minRecombinant(*it2);
    }
  }
  cout << endl;

  return (0);
}

int pedhap::genotypingError(string outfile)
{
  ofstream outf(outfile.c_str());
  cout << "Writing probable genotype errors to " << outfile << endl;
  outf << "ID\tfather\tmother\tRSID1\tPOS\tPROB_ERROR\n";
  for (vector<set<string>>::iterator it1 = ped->pedigrees.begin(); it1 != ped->pedigrees.end(); it1++)
  {

    vector<string> idorder;
    ped->orderSamples(*it1, idorder);

    for (vector<string>::iterator it2 = idorder.begin(); it2 != idorder.end(); it2++)
    {
      string dad = ped->sampleinfo[*it2].dad;
      string mum = ped->sampleinfo[*it2].mum;
      cout << "Detecting errors for: " << *it2 << " - " << dad << " - " << mum << endl;

      if (dad.compare("0") != 0 && mum.compare("0") != 0)
      {

        vector<bool> *c = haps->getHap(*it2);
        vector<bool> *d = haps->getHap(dad);
        vector<bool> *m = haps->getHap(mum);
        trio->setHaps(d, m, c);
        trio->EM(NITERATION);
        for (int l = 0; l < nsnp; l++)
        {
          float p = trio->genError[l];
          if (p > 0.1)
          {
            outf << *it2 << "\t" << dad << "\t" << mum;
            outf << "\t" << haps->rsid1[l] << "\t" << haps->positions[l];
            outf << "\t" << trio->genError[l] << endl;
          }
        }
      }
      else if (dad.compare("0") != 0)
      {

        vector<bool> *c = haps->getHap(*it2);
        vector<bool> *d = haps->getHap(dad);
        duo->setHaps(d, c, "1");
        duo->EM(NITERATION);
        for (int l = 0; l < nsnp; l++)
        {
          float p = duo->genError[l];
          if (p > 0.1)
          {
            outf << *it2 << "\t" << dad << "\t" << mum;
            outf << "\t" << haps->rsid1[l] << "\t" << haps->positions[l];
            outf << "\t" << duo->genError[l] << endl;
          }
        }
      }
      else if (mum.compare("0") != 0)
      {

        vector<bool> *c = haps->getHap(*it2);
        vector<bool> *m = haps->getHap(mum);
        duo->setHaps(m, c, "2");
        duo->EM(NITERATION);
        for (int l = 0; l < nsnp; l++)
        {
          float p = duo->genError[l];
          if (p > 0.1)
          {
            outf << *it2 << "\t" << dad << "\t" << mum;
            outf << "\t" << haps->rsid1[l] << "\t" << haps->positions[l];
            outf << "\t" << duo->genError[l] << endl;
          }
        }
      }
    }
  }

  return 0;
}

int pedhap::recombinationMap(string outfile)
{
  ofstream outf(outfile.c_str());

  outf << "PARENT CHILD";
  for (int l = 0; l < nsnp; l++)
    outf << " " << haps->positions[l];
  outf << endl;

  cout << "Calculating recombination map and writing to " << outfile << endl;
  for (vector<set<string>>::iterator it1 = ped->pedigrees.begin(); it1 != ped->pedigrees.end(); it1++)
  {

    vector<string> idorder;
    ped->orderSamples(*it1, idorder);

    for (vector<string>::iterator it2 = idorder.begin(); it2 != idorder.end(); it2++)
    {
      string dad = ped->sampleinfo[*it2].dad;
      string mum = ped->sampleinfo[*it2].mum;

      if (dad.compare("0") != 0 || mum.compare("0") != 0)
      {

        cout << "Mapping recombination for: " << *it2 << " - " << dad << " - " << mum << endl;
        if (dad.compare("0") != 0 && mum.compare("0") != 0)
        {
          vector<bool> *c = haps->getHap(*it2);
          vector<bool> *d = haps->getHap(dad);
          vector<bool> *m = haps->getHap(mum);
          trio->setHaps(d, m, c);
          trio->estimateRecombination();
          outf << *it2 << " " << dad << "\t";
          for (int l = 0; l < nsnp; l++)
            outf << trio->recombinationPat[l] << " ";
          outf << endl
               << *it2 << " " << mum << "\t";
          for (int l = 0; l < nsnp; l++)
            outf << trio->recombinationMat[l] << " ";
          outf << endl;
        }
        else if (dad.compare("0") != 0)
        {

          vector<bool> *c = haps->getHap(*it2);
          vector<bool> *d = haps->getHap(dad);
          duo->setHaps(d, c, "1");
          duo->estimateRecombination();
          outf << *it2 << " " << dad << "\t";
          for (int l = 0; l < nsnp; l++)
            outf << duo->recombinationMap[l] << " ";
          outf << endl;
        }
        else if (mum.compare("0") != 0)
        {
          vector<bool> *c = haps->getHap(*it2);
          vector<bool> *m = haps->getHap(mum);
          duo->setHaps(m, c, "2");
          duo->estimateRecombination();
          outf << *it2 << " " << mum << "\t";
          for (int l = 0; l < nsnp; l++)
            outf << duo->recombinationMap[l] << " ";
          outf << endl;
        }
      }
    }
  }
  return 0;
}


