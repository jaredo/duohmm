#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "hapmodule.h"

int main(int argc, char **argv)
{
    assert(argc==4);
    Haplotypes haps(argv[1]);
    haps.setMissing(argv[2]);
    int count =         haps.imputeAll();
    cerr<<  "Re-imputed "<< count <<" genotypes"<<endl;
    haps.writeHaps(argv[3]);
    return(0);
}