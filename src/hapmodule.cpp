#include "hapmodule.h"

#define DEBUG 0
Haplotypes::Haplotypes(const Haplotypes& h) {
  cout << "deep copying"<<endl;
  nhap = h.nhap;
  nsnp = h.nsnp;

  H = newMatrix<unsigned char>(nhap,nsnp);
  for(int i=0;i<nhap;i++) memcpy(H[i],h.H[i],nsnp);
}


Haplotypes::Haplotypes(string filename) {
  input_file = filename;
  ifstream inf1((filename+".sample").c_str());
  string line;
  int pos;
  char c;
  int count = 0;
  //int *positions;
  ids.clear();
  string id;
  while( getline(inf1,line) ) {
    if(count>1) {
      istringstream iss(line);
      iss >> id;iss >> id;
      ids.push_back(id);
      idlook[id] = count-2;
      //      cout << count << ":\t"<< ids[count-2] << endl;
    }
    count++;
  }
  cout << count-2 << " individuals in " << filename <<".sample"<< endl;
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
  cout << count << " markers in " << filename <<".haps" << endl;
  nsnp = count;
  positions.resize(nsnp);
  rsid1.resize(nsnp);
  rsid2.resize(nsnp);
  ref.resize(nsnp);
  alt.resize(nsnp);

  blah.seekg(0);
  blah.clear();
  inf2.pop();
  inf2.push(blah);

  H = newMatrix<unsigned char>(nhap,nsnp);
  string tmp;
  for(int i=0;i<nsnp;i++) {
    if(!inf2) {
      cerr << "Problem reading "<<filename<<".haps" <<endl;
      exit(1);
    }
    inf2 >> tmp;
    rsid1[i]=tmp;
    inf2 >> tmp;
    rsid2[i]=tmp;
    inf2 >> pos;
    positions[i] = pos;
    inf2 >> tmp;
    ref[i]=tmp;
    inf2 >> tmp;
    alt[i]=tmp;
    //    cout << rsid1[i] << " "<< rsid2[i] << " " << positions[i] << " " << ref[i] << " " << alt[i] << endl;
    //for(int j=0;j<5;j++) inf2.ignore(1000,' ');
    for(int j=0;j<nhap;j++) {
      unsigned   int hapval;
      //      inf2 >> c;
      //      H[j][i] = atoi(&c);
      inf2 >> hapval;
      H[j][i] = hapval;
      if(!(H[j][i]==0 || H[j][i]==1)) {
	cerr << "INVALID VALUE IN HAPLOTYPE FILE ON LINE " << i+1 << " COLUMN " << j+6 <<  "\n" << hapval << " " <<(unsigned int) H[j][i] << "\nExiting..."<<endl;
	exit(1);
      }
      if(DEBUG>0)      if(j<4 && i<20) cout << (int)H[j][i] << " ";
    }
    if(DEBUG>0)    if(i<20) cout << endl;
  }

  cout << nhap << " haplotypes over " << nsnp << " markers." << endl;
}

int Haplotypes::writeHaps(string fname) {
  string outfname = fname + ".haps";
  if(fileexists(outfname)) {
    cerr << "ERROR: " << outfname<<" exists.  Will not overwrite.\nExiting..."<<endl;
    exit(1);
  }

  cout << "Writing pedigree corrected haplotypes to " << outfname << endl;
  ofstream outf(outfname.c_str());
  for(int l=0;l<nsnp;l++) {
    outf<<rsid1[l]<<" "<<rsid2[l]<<" "<<positions[l]<<" "<<ref[l]<<" "<<alt[l];
    for(int i=0;i<nhap/2;i++) {
      outf << " "<< (int)H[i*2][l];
      outf << " "<< (int)H[i*2+1][l];
    }
    outf<<endl;
  }

  string infname=input_file+".sample";
  outfname = fname+".sample";
  ifstream  src(infname.c_str(), std::ios::binary);
  if(!src) {
    cerr << "Probleam reading "<<infname <<endl;
    exit(1);
  }
  ofstream  dst(outfname.c_str(),   std::ios::binary);
  dst << src.rdbuf();

  return(0);
}


Haplotypes::~Haplotypes() {
  delMatrix<unsigned char>(H,nhap,nsnp);
}

unsigned char **Haplotypes::getHap(string id){
  assert(idlook.count(id));
  return(&H[idlook[id]*2]);
}

