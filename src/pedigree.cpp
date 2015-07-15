#include "pedigree.h"

#define DEBUG 0

pedigree::pedigree(string fname,vector<string> & ids,char type) {
  if(type=='f')
    fromFam(fname,ids);
  else if(type=='s')
    fromSample(fname,ids);
}

int pedigree::fromSample(string fname,vector<string> & ids) {
  io::filtering_istream inf2;
  ifstream blah(fname.c_str(),ios_base::in);
  inf2.push(blah);
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
  delete header;
  return 0;
}



int pedigree::fromFam(string fname,vector<string> & ids) {
  io::filtering_istream inf2;
  ifstream blah(fname.c_str(),ios_base::in);
  inf2.push(blah);
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

  cout << pedigrees.size() << " distinct pedigrees found."<<endl;
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

  cout << "SIZE\tFREQUENCY"<<endl;
  for (map<int,int>::iterator it2=freq.begin(); it2!=freq.end(); it2++)
    cout << it2->first << "\t" << it2->second << endl;


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

