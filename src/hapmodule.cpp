#include "hapmodule.h"

Haplotypes::Haplotypes(const Haplotypes& h) {
	cout << "deep copying"<<endl;
	nhap = h.nhap;
	nsnp = h.nsnp;

	H = newMatrix<unsigned char>(nhap,nsnp);
	for(int i=0;i<nhap;i++) memcpy(H[i],h.H[i],nsnp);

}


Haplotypes::Haplotypes(string filename) {
	ifstream inf1((filename+".sample").c_str());
	string line,rsid1,rsid2,pos,ref,alt;
	char c;
	int count = 0;
	//int *positions;

	while ( getline(inf1,line) )
		count ++;
	cout << count << " lines in " << filename <<".sample"<< endl;
	nhap = (count-2) * 2;
	allsamples.resize(nhap);
	for(int i=0;i<nhap;i++) allsamples[i] = i;
	inf1.close();


	io::filtering_istream inf2;
	//ifstream blah((filename+".haps").c_str(),ios_base::in | ios_base::binary);
	ifstream blah((filename+".haps").c_str(),ios_base::in);
		//inf2.push(io::gzip_decompressor());
	inf2.push(blah);
	count = 0;
	while ( getline(inf2,line) ) {
		count ++;
	}
	cout << count << " lines in " << filename <<".haps" << endl;
	nsnp = count;

	blah.seekg(0);
	blah.clear();
	inf2.pop();
	inf2.push(blah);

	H = newMatrix<unsigned char>(nhap,nsnp);

	for(int i=0;i<nsnp;i++) {

		if(!inf2) {
			cerr << "Problem reading "<<filename<<".haps" <<endl;
			exit(1);
		}
		 inf2 >> rsid1;
		 inf2 >> rsid2;
		 inf2 >> pos;
		 inf2 >> ref;
		 inf2 >> alt;
//		 cout << rsid1 << " "<< rsid2 << " " << pos << " " << ref << " " << alt << endl;
		for(int j=0;j<5;j++) inf2.ignore(1000,' ');
		for(int j=0;j<nhap;j++) {
			inf2 >> c;
			H[j][i] = (unsigned char)atoi(&c);
		}
	}

	cout << nhap << " haplotypes over " << nsnp << " markers." << endl;
}

Haplotypes::~Haplotypes() {
	delMatrix<unsigned char>(H,nhap,nsnp);
}
