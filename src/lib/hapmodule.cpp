//$Id$

#include "hapmodule.h"

#define DEBUG 0

Haplotypes::~Haplotypes()
{
  delete[] H;
}

bool Haplotypes::isMissing(string &id, int index)
{
  assert(!_missing.empty());
  assert(idlook.count(id));
  return _missing[idlook[id]][index];
}

vector<bool> *Haplotypes::getHap(string id)
{
  assert(idlook.count(id));
  return (&H[idlook[id] * 2]);
}

int Haplotypes::setMissing(const string &fname)
{
  ifstream inf1(fname.c_str());
  if (!inf1)
  {
    cerr << "Problem reading " << fname << endl;
    exit(1);
  }
  string line;
  string sampleid, chrom, rsid;
  int pos, count = 0;
  _missing.assign(nhap / 2, vector<bool>(nsnp, false));
  map<int,int> poslook;
  for(size_t i=0;i<positions.size();i++)
    poslook[positions[i]]=i;

  while (getline(inf1, line))
  {
    istringstream ss(line);
    ss >> sampleid;
    ss >> chrom;
    ss >> rsid;
    ss >> pos;
    _missing[idlook[sampleid]][poslook[pos]] = true;
    count++;
  }
  cerr << "Read " << count << " missing genotypes from " << fname << std::endl;
  return (count);
}

int Haplotypes::imputeAll()
{
  int count = 0;
  for(auto s=ids.begin();s!=ids.end();s++)
  {
    count+=imputeSample(*s);
  }    
  return count;
}

int Haplotypes::imputeSample(string const & id)
{
  int sample_index = idlook[id];
  int count=0;
  for(size_t i=0;i<nsnp;i++)
  {
      if(_missing[sample_index][i])
      {
          float i0= impute(sample_index*2,i);
          float i1= impute(sample_index*2+1,i);
          H[sample_index*2][i] = round(i0);
          H[sample_index*2+1][i] = round(i1);
         // fprintf(stdout,"%s\t%d\t%d/%d\t%.1f/%.1f\t%d\t%d\n",id.c_str(),positions[i],g0,g1,i0,i1,g0+g1,(int)round(i0+i1));
          count++;
      }
  }    
  return count;
}

float Haplotypes::impute(int hap_index,int pos_index)
{
  Lmatch.assign(nhap,0.);
  Rmatch.assign(nhap,0.);
  float num=0,den=0;
  for(size_t i=0;i<nhap;i++)
  {
    if(_missing[i/2][pos_index] || i==hap_index)
    {
      Lmatch[i]=0.;
      Rmatch[i]=0.;
//      std::cerr << "skipped "<<ids[i/2]<<endl;
    }
    else
    {
      if(pos_index==nsnp-1)
      {
        Rmatch[i]=1.;
      }
      else{
        int cursor = pos_index+1;
        while(cursor<nsnp && H[i][cursor]==H[hap_index][cursor])    cursor++;
        if(cursor==nsnp)
          Rmatch[i] = positions[nsnp-1] - positions[pos_index];
        else        
          Rmatch[i] = (positions[cursor]+positions[cursor-1])/2 - positions[pos_index];

        Rmatch[i] /= 1e3;
      }
      if(pos_index==0)
      {
        Lmatch[i]=1.;
      }
      else {
        int cursor = pos_index-1;
        while(cursor>=0 && H[i][cursor]==H[hap_index][cursor]) cursor--;
        if(cursor==-1)
          Lmatch[i] = positions[pos_index] - positions[0];
        else
          Lmatch[i] = positions[pos_index] - (positions[cursor]+positions[cursor+1])/2;
        Lmatch[i] /= 1e3;
      }
      num += H[i][pos_index]*Lmatch[i]*Rmatch[i];
      den += Lmatch[i]*Rmatch[i];
    }
  }
  if(den>0)
  {
    return num/den;
  }    
  else
  {
    cerr << "WARNING: numeric issue at " <<hap_index <<" "<<positions[pos_index]<<endl;
    return (float)H[hap_index][pos_index];
  }    
}

Haplotypes::Haplotypes(string filename)
{
  input_file = filename;
  ifstream inf1((filename + ".sample").c_str());
  if (!inf1)
  {
    cerr << "Problem reading " << filename + ".sample" << endl;
    exit(1);
  }

  string line;
  int pos;
  int count = 0;
  //int *positions;
  ids.clear();
  string id;
  while (getline(inf1, line))
  {
    if (count > 1)
    {
      istringstream iss(line);
      iss >> id;
      iss >> id;
      ids.push_back(id);
      idlook[id] = count - 2;
      //      cout << count << ":\t"<< ids[count-2] << endl;
    }
    count++;
  }
  cout << count - 2 << " individuals in " << filename << ".sample" << endl;
  nhap = (count - 2) * 2;
  allsamples.resize(nhap);
  for (int i = 0; i < nhap; i++)
    allsamples[i] = i;
  inf1.close();
  //ifstream blah((filename+".haps").c_str(),ios_base::in | ios_base::binary);
  ifstream inf2((filename + ".haps").c_str());
  if (!inf2)
  {
    cerr << "Problem reading " << filename + ".haps" << endl;
    exit(1);
  }

  count = 0;
  while (getline(inf2, line))
  {
    count++;
  }
  cout << count << " markers in " << filename << ".haps" << endl;
  nsnp = count;
  positions.resize(nsnp);
  rsid1.resize(nsnp);
  rsid2.resize(nsnp);
  ref.resize(nsnp);
  alt.resize(nsnp);

  inf2.close();
  inf2.open((filename + ".haps").c_str());

  H = new vector<bool>[nhap];
  for (int j = 0; j < nhap; j++)
    H[j].assign(nsnp, 0);
  string tmp;
  for (int i = 0; i < nsnp; i++)
  {
    if (!inf2)
    {
      cerr << "Problem reading " << filename << ".haps at line " << i << endl;
      exit(1);
    }
    inf2 >> tmp;
    rsid1[i] = tmp;
    inf2 >> tmp;
    rsid2[i] = tmp;
    inf2 >> pos;
    positions[i] = pos;
    inf2 >> tmp;
    ref[i] = tmp;
    inf2 >> tmp;
    alt[i] = tmp;
    //    cout << rsid1[i] << " "<< rsid2[i] << " " << positions[i] << " " << ref[i] << " " << alt[i] << endl;
    //for(int j=0;j<5;j++) inf2.ignore(1000,' ');
    for (int j = 0; j < nhap; j++)
    {
      int hapval;
      //      inf2 >> c;
      //      H[j][i] = atoi(&c);
      inf2 >> hapval;
      H[j][i] = hapval;
      if (!(hapval == 0 || hapval == 1))
      {
        cerr << "INVALID VALUE IN HAPLOTYPE FILE ON LINE " << i + 1 << " COLUMN " << j + 6 << "\n"
             << hapval << " " << (unsigned int)H[j][i] << "\nExiting..." << endl;
        exit(1);
      }
      if (DEBUG > 0)
        if (j < 4 && i < 20)
          cout << (int)H[j][i] << " ";
    }
    if (DEBUG > 0)
      if (i < 20)
        cout << endl;
  }

  cout << nhap << " haplotypes over " << nsnp << " markers." << endl;
}

int Haplotypes::writeHaps(string fname)
{
  string outfname = fname + ".haps";
  if (fileexists(outfname))
  {
    cerr << "ERROR: " << outfname << " exists.  Will not overwrite.\nExiting..." << endl;
    exit(1);
  }

  cout << "Writing pedigree corrected haplotypes to " << outfname << endl;
  ofstream outf(outfname.c_str());
  for (int l = 0; l < nsnp; l++)
  {
    outf << rsid1[l] << " " << rsid2[l] << " " << positions[l] << " " << ref[l] << " " << alt[l];
    for (int i = 0; i < nhap / 2; i++)
    {
      outf << " " << (int)H[i * 2][l];
      outf << " " << (int)H[i * 2 + 1][l];
    }
    outf << endl;
  }

  string infname = input_file + ".sample";
  outfname = fname + ".sample";
  ifstream src(infname.c_str(), std::ios::binary);
  if (!src)
  {
    cerr << "Problem reading " << infname << endl;
    exit(1);
  }
  ofstream dst(outfname.c_str(), std::ios::binary);
  dst << src.rdbuf();

  return (0);
}

#ifdef SHAPEIT
//these routines deals specifically with converting to/from SHAPEIT2 data structures

Haplotypes::Haplotypes(filter_writer &F, genhap_set &GH)
{ //builds haps from a S2 genhaptset
  _F = &F;
  _GH = &GH;
  //haplotype_writer::haplotype_writer(filter_writer & _F, genhap_set & _GH, string header1, string header2, bool founder) : F(_F), GH(_GH) {

  //INTERNAL DATA

  nhap = GH.vecG.size() * 2;
  nsnp = GH.mapG->size();
  if (DEBUG > 0)
    cout << nhap << " haps and " << nsnp << " snps" << endl;
  //this bit orders the data...i think
  for (int i = 0; i < GH.vecG.size(); i++)
  {
    genotype_graph *g = GH.vecG[i];
    int type = g->type;
    assert(!(type == DUO_F || type == TRIO_F) && !(type == DUO_M || type == TRIO_M));
    if (F.checkInd(g->name))
      orderI.push_back(haplotype_index(g->index, i, UNR));
  }

  assert(orderI.size() == nhap / 2);
  //now we label the samples
  ids.resize(nhap / 2);
  for (int i = 0; i < orderI.size(); i++)
  {
    genotype_graph *g = GH.vecG[orderI[i].indexG];
    int type = orderI[i].type;
    string name;
    if (type == DUO_F || type == TRIO_F)
    {
      name = g->fname;
    }
    else if (type == DUO_M || type == TRIO_M)
    {
      name = g->mname;
    }
    else
    {
      name = g->name;
    }
    ids[i] = name;
    idlook[name] = i;
  }

  H = new vector<bool>[nhap];
  for (int i = 0; i < nhap; i++)
    H[i].assign(nsnp, 0);
  positions.resize(nsnp);
  ref.assign(nsnp, "A");
  alt.assign(nsnp, "A");
  rsid1.resize(nsnp, "-");
  rsid2.resize(nsnp, "-");
  _missing.assign(nhap / 2, vector<bool>(nsnp, false));
  //now we convert the haplotypes.
  for (int l = 0; l < GH.mapG->size(); l++)
  {
    snp *s = GH.mapG->vec_pos[l];
    if (F.checkSnp(s->bp))
    {
      positions[l] = s->bp;
      for (int i = 0; i < orderI.size(); i++)
      {
        int type = orderI[i].type;
        unsigned int a1, a2;
        assert(!(type == DUO_F || type == TRIO_F) && !(type == DUO_M || type == TRIO_M));
        _missing[i][l] = GH.vecG[orderI[i].indexG]->Miss[l];
        a1 = GH.vecH[GH.G2H[orderI[i].indexG][0]][l];
        a2 = GH.vecH[GH.G2H[orderI[i].indexG][1]][l];

        H[i * 2][l] = a1;
        H[i * 2 + 1][l] = a2;
      }
    }
  }
  allsamples.resize(nhap / 2);
  for (int i = 0; i < nhap / 2; i++)
    allsamples[i] = i;
  if (DEBUG > 0)
    cout << "finished reading haps " << endl;
}

int Haplotypes::getSHAPEIT2(filter_writer &F, genhap_set &GH)
{ //puts haplotypes back  into GH for S2
  for (int l = 0; l < GH.mapG->size(); l++)
  {
    snp *s = GH.mapG->vec_pos[l];
    if (F.checkSnp(s->bp))
    {
      positions[l] = s->bp;
      for (int i = 0; i < orderI.size(); i++)
      {
        int type = orderI[i].type;
        assert(!(type == DUO_F || type == TRIO_F) && !(type == DUO_M || type == TRIO_M));
        GH.vecH[GH.G2H[orderI[i].indexG][0]].set(l, H[i * 2][l]);
        GH.vecH[GH.G2H[orderI[i].indexG][1]].set(l, H[i * 2 + 1][l]);
      }
    }
  }
  return (0);
}


vector<bool> *Haplotypes::getMissing(string &id)
{
  if (idlook.count(id))
  {
    return (&_missing[idlook[id]]);
  }
  else
  {
    return nullptr;
  }
}

#endif
