#include "pedigree.h"

#define DEBUG 0

pedigree::pedigree(string fname,vector<string> & ids,char type) {
  if(type=='f')
    fromFam(fname,ids);
  else if(type=='s')
    fromSample(fname,ids);
}

int pedigree::fromSample(string fname,vector<string> & ids) {
  ifstream inf2(fname.c_str());
  set<string> idcheck(ids.begin(),ids.end());
  string fid,id,dad,mum,sex;
  int count = 0;
  char *header = new char[10000];
  inf2.getline(header,10000);
  istringstream iss(header);
  string field;
  int dadidx=-1,mumidx=-1,sexidx=-1;

  count=0;
  while(iss) {
    iss >> field;
    if(field=="father"||field=="ID_father")
      dadidx=count;
    if(field=="mother"||field=="ID_mother")
      mumidx=count;
    if(field.compare("sex")==0)
      sexidx=count;
    count++;
  }
  if(DEBUG>0) {
    cout << "dadidx = "<< dadidx<<endl;
    cout << "mumidx = "<< mumidx<<endl;
    cout << "sexidx = "<< sexidx<<endl;
  }
  //  inf2.ignore('\n',10000);

  if(dadidx==-1) {
    cerr<<"ERROR: father column was not in "<<fname<<endl;
    exit(1);
  }
  if(mumidx==-1) {
    cerr<<"ERROR: mother column was not in "<<fname<<endl;
    exit(1);
  }
  if(sexidx==-1) {
    cerr<<"ERROR: sex column was not in "<<fname<<endl;
    exit(1);
  }

  vector<string> line(count);
  int ncol = count-1;
  cout << ncol << " columns in " << fname << endl;
  //  cout << header << endl;
  inf2.getline(header,10000);
  //  cout << header << endl;


  while(inf2) {
    for(int i=0;i<ncol;i++) 
      inf2 >> line[i];

    individual tmp;
    tmp.fid = line[0];
    tmp.id = line[1];
    tmp.dad = line[dadidx];
    tmp.mum = line[mumidx];
    tmp.sex = line[sexidx];
    if(DEBUG>0) cout << tmp.fid << "\t" << tmp.id << "\t" << tmp.dad << "\t" << tmp.mum << "\t" << tmp.sex << endl;

    if(tmp.id.compare("0")==0) {//checks there are no individuals with "0" as id
      cerr << "ERROR: Individual has ID=0\nExiting..."<<endl;
      exit(1);
    }

    if(idcheck.count(tmp.dad)==0) tmp.dad = "0";
    if(idcheck.count(tmp.mum)==0) tmp.mum = "0";

    tmp.idx = count;
    tmp.fidx = -1;
    inf2.ignore(10000,'\n');

    if(inf2) {
      count++;
      if(idcheck.count(tmp.id)>0) sampleinfo[tmp.id] = tmp;
    }
    
  }
  cout << count << " samples in " << fname << endl;
  buildPeds();
  delete[] header;
  return 0;
}



int pedigree::fromFam(string fname,vector<string> & ids) {
  ifstream inf2(fname.c_str());
  set<string> idcheck(ids.begin(),ids.end());
  string fid,id,dad,mum,sex;
  int count = 0;

  while(inf2) {
    individual tmp;
    inf2 >> tmp.fid;
    inf2 >> tmp.id;
    inf2 >> tmp.dad;
    inf2 >> tmp.mum;
    inf2 >> tmp.sex;

    if(tmp.id.compare("0")==0) {//checks there are no individuals with "0" as id
      cerr << "ERROR: Individual has ID=0\nExiting..."<<endl;
      exit(1);
    }

    if(idcheck.count(tmp.dad)==0) tmp.dad = "0";
    if(idcheck.count(tmp.mum)==0) tmp.mum = "0";

    tmp.idx = count;
    tmp.fidx = -1;
    inf2.ignore(1000,'\n');
    if(inf2) {
      count++;
      if(idcheck.count(tmp.id)>0) sampleinfo[tmp.id] = tmp;
      //			cout << count << " " << tmp.id << endl;
    }
  }
  cout << count << " samples in " << fname << endl;
  buildPeds();
  return 0;
}

int buildFam(string id, map<string,individual> & sampleinfo, set<string> & fam) {
  //  cout << "ID " << id << endl;
  if(!fam.count(id) && sampleinfo.count(id)) {
    fam.insert(id);
    if(sampleinfo.count(sampleinfo[id].dad))
      buildFam(sampleinfo[id].dad,sampleinfo,fam);

    if(sampleinfo.count(sampleinfo[id].mum))
      buildFam(sampleinfo[id].mum,sampleinfo,fam);

    for (map<string,individual>::iterator it=sampleinfo.begin(); it!=sampleinfo.end(); it++) {
      if(it->second.dad.compare(id)==0) {
	if( sampleinfo[id].sex.compare("1")!=0 )
	  cerr << "WARNING: sex of sample "<<id<<" does not match pedigree information.\n";
	buildFam(it->second.id,sampleinfo,fam);
	sampleinfo[id].kids.insert(it->second.id);
      }
      if(it->second.mum.compare(id)==0) {
	if(sampleinfo[id].sex.compare("2")!=0)
	  cerr << "WARNING: sex of sample "<<id<<" does not match pedigree information.\n";
	buildFam(it->second.id,sampleinfo,fam);
	sampleinfo[id].kids.insert(it->second.id);
      }
    }
  }
  return 0;
}

int pedigree::buildPeds() {

  for (map<string,individual>::iterator it1=sampleinfo.begin(); it1!=sampleinfo.end(); it1++) {
    bool hasfamily = false;

    for (vector<set<string> >::iterator it2=pedigrees.begin(); it2!=pedigrees.end(); it2++)
      if(it2->count(it1->first))
	hasfamily=true;

    if(!hasfamily) {
      set<string> newfam;
      buildFam(it1->first,sampleinfo,newfam);
      pedigrees.push_back(newfam);
    }
  }

  map<int,int> freq;
  for (vector<set<string> >::iterator it2=pedigrees.begin(); it2!=pedigrees.end(); it2++) {
    int n = it2->size();
    if(!freq.count(n)) freq[n] = 1;
    else freq[n]++;

    if(DEBUG>1) {
      if(n>0) {
	for(set<string>::iterator it3=it2->begin(); it3!=it2->end(); it3++)
	  cout << *it3 << " ";
	cout << endl;
      }
    }

  }

#ifndef SHAPEIT
  cout << pedigrees.size() << " distinct pedigrees found."<<endl;
  cout << "SIZE\tFREQUENCY"<<endl;
  for (map<int,int>::iterator it2=freq.begin(); it2!=freq.end(); it2++)
    cout << it2->first << "\t" << it2->second << endl;
#endif

  return(0);
}


int pedigree::orderSamples(set<string> & ids,vector<string> & ordered_ids) {

  if(ids.size()==1) {
    ordered_ids.push_back(*ids.begin());
    ids.erase(*ids.begin());
  }

  if(ids.size()>1) {
    set<string> newids;
    for(set<string>::iterator it1=ids.begin(); it1!=ids.end(); it1++) {
      if( (sampleinfo[*it1].dad.compare("0")==0||!ids.count(sampleinfo[*it1].dad)) && (sampleinfo[*it1].mum.compare("0")==0||!ids.count(sampleinfo[*it1].mum)))
	ordered_ids.push_back(*it1);
      else
	newids.insert(*it1);
    }
    orderSamples(newids,ordered_ids);
  }

  return(0);
}

#ifdef SHAPEIT

pedigree::pedigree(filter_writer & F, genhap_set & GH, string header1, string header2,vector<string> & ids) {
  if(DEBUG>0) cout << header1<<endl<<header2<<endl;
  istringstream iss(header1);
  string field;
  int count = 0;
  set<string> idcheck(ids.begin(),ids.end());
  int dadidx=-1,mumidx=-1,sexidx=-1;
  while(iss) {
    iss >> field;
    if(!iss)
      break;
    if(field.compare("father")==0||field.compare("ID_father")==0)
      dadidx=count;
    if(field.compare("mother")==0||field.compare("ID_mother")==0)
      mumidx=count;
    if(field.compare("sex")==0||field.compare("gender")==0)      
      sexidx=count;
    count++;
  }
  dadidx--;
  mumidx--;
  sexidx--;

  if(DEBUG>0) cout << dadidx << " " << mumidx << " " << sexidx << endl;
  count = 0;
  bool sex_warn=false;
  for (int i = 0 ; i < GH.vecG.size(); i ++) {
    genotype_graph * g = GH.vecG[i];
    int type = g->type;
    assert(!(type == DUO_F || type == TRIO_F)  && !(type == DUO_M || type == TRIO_M));
    /*
    if (type == DUO_F || type == TRIO_F) {
      if(DEBUG>0) cout  << g->fproperties[0] << " " << g->fname;
      for (int p = 1 ; p < g->fproperties.size() ; p ++) if(DEBUG>0) cout  << " " << g->fproperties[p];
      if(DEBUG>0) cout  << endl;
    } else if (type == DUO_M || type == TRIO_M) {
      if(DEBUG>0) cout  << g->mproperties[0] << " " << g->mname;
      for (int p = 1 ; p < g->mproperties.size() ; p ++) if(DEBUG>0) cout  << " " << g->mproperties[p];
      if(DEBUG>0) cout  << endl;
    } else {
      if(DEBUG>0) cout  << g->properties[0] << " " << g->name;
      for (int p = 1 ; p < g->properties.size() ; p ++) if(DEBUG>0) cout  << " " << g->properties[p];
      if(DEBUG>0) cout  << endl;
    }
    */
    if(!(dadidx < g->properties.size() && dadidx>=0)) {
      cerr <<"ERROR: duohmm module could not find father id column"<<endl;
      exit(1);
    }
    if(!(mumidx < g->properties.size() && mumidx>=0)) {
      cerr <<"ERROR: duohmm module could not find mother id column"<<endl;
      exit(1);      
    }
    if(!sex_warn && !(sexidx < (int)g->properties.size() && sexidx>=0))  {
      cerr << "WARNING: sample sex was not provided"<<endl;
      sex_warn=true;
    }
    
    individual tmp;
    tmp.fid = g->properties[0];
    tmp.id =  g->name;
    tmp.dad = g->properties[dadidx];
    tmp.mum = g->properties[mumidx];
    if(sexidx>=0)
      tmp.sex = g->properties[sexidx];
    else
      tmp.sex="unknown";

    if(DEBUG>0)
      cout << tmp.fid << " " << tmp.id << " " << tmp.dad << " " <<tmp.mum << " " <<tmp.sex<<endl;

    assert(tmp.id.compare("0")!=0);

    if(idcheck.count(tmp.dad)==0) tmp.dad = "0";
    if(idcheck.count(tmp.mum)==0) tmp.mum = "0";

    tmp.idx = count;
    tmp.fidx = -1;
    count++;
    if(idcheck.count(tmp.id)>0) sampleinfo[tmp.id] = tmp;
  }
  if(DEBUG>0) cout << count << " samples" << endl;
  buildPeds();
}

#endif
