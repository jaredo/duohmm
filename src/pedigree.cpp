#include "pedigree.h"


pedigree::pedigree(string fname){
	io::filtering_istream inf2;
	ifstream blah(fname.c_str(),ios_base::in);
	inf2.push(blah);

	string fid,id,dad,mum,sex;
	int count = 0;

	while(inf2) {
		individual tmp;
		inf2 >> tmp.fid;
		inf2 >> tmp.id;
		inf2 >> tmp.dad;
		inf2 >> tmp.mum;
		inf2 >> tmp.sex;

		assert((tmp.id.compare("0")!=0));//checks there are no individuals with "0" as id

		tmp.idx = count;
		inf2.ignore(1000,'\n');
		if(inf2) {
			count++;
			sampleinfo[tmp.id] = tmp;
			//			cout << count << " " << tmp.id << endl;
		}
	}
	cout << count << " pedigrees in " << fname << endl;
	buildPeds();
}

int buildFam(string id, map<string,individual> & sampleinfo, set<string> & fam) {
	if(!fam.count(id) && sampleinfo.count(id)) {
		fam.insert(id);
		if(sampleinfo.count(sampleinfo[id].dad))
			buildFam(sampleinfo[id].dad,sampleinfo,fam);

		if(sampleinfo.count(sampleinfo[id].dad))
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
		if(n>1) {
			for(set<string>::iterator it3=it2->begin(); it3!=it2->end(); it3++)
				cout << *it3 << " ";
			cout << endl;
		}
	}

	cout << "SIZE\tFREQUENCY"<<endl;
	for (map<int,int>::iterator it2=freq.begin(); it2!=freq.end(); it2++)
		cout << it2->first << "\t" << it2->second << endl;


	return 0;
}


