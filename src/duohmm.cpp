//============================================================================
// Name        : duohmm.cpp
// Author      : Jared O'Connell
// Copyright   : Jared O'Connell (2014)
// Description : post-hoc pedigree integration for shapeit2 haplotypes.
//============================================================================

#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include "hapmodule.h"
#include "pedigree.h"
#include "hmm.h"
#include "pedhap.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

int main(int argc, char **argv)
{
  try
  {
    string haps, fam, gm, corrected_out, generr_out, rec_out;
    int niteration;
    po::options_description desc("Allowed options");
    desc.add_options()
        //      ("help", "produce help message")
        ("haps,H", po::value<string>(&haps), "SHAPEIT2 haplotype file")
        //      ("fam,F", po::value<string>(&fam), "pedigree information in PLINK .fam format")
        ("input-gen,M", po::value<string>(&gm)->default_value(""), "genetic map file in SHAPEIT2/Impute format")("output-hap,O", po::value<string>(&corrected_out), "output pedigree corrected haplotypes this file")("output-generr,G", po::value<string>(&generr_out), "output possible genotyping errors to this file")("output-rec,R", po::value<string>(&rec_out), "output recombination map for pedigree to this file")("niteration", po::value<int>(&niteration)->default_value(100), "maximum number of EM iterations for parameter estimation");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    cout << endl;
    cout << "DuoHMM: a post-hoc method for integrating extended pedigree information with SHAPEIT2 haplotypes\n"
         << endl;
    cout << "\t* Authors: Jared O'Connell, Olivier Delaneau and Jonathan Marchini" << endl;
    cout << "\t* Webpage: http://www.stats.ox.ac.uk/~marchini/software/gwas/gwas.html" << endl;
    cout << "\t* Version: " << VERSION << endl;
    cout << "\t* Contact: jaredoconnell@gmail.com" << endl;
    cout << endl;

    if (vm.count("help") || argc == 1)
    {
      cout << desc << "\n";
      return 0;
    }

    if (vm.count("haps"))
      cout << "Haplotypes file:\t" << vm["haps"].as<string>() << "\n";
    else
    {
      cout << "No haplotypes file provided.\nExiting..." << endl;
      return (1);
    }
    if (!(vm.count("output-hap") || vm.count("output-generr") || vm.count("output-rec")))
    {
      cout << "You must specify one of: --output-hap --output-generr --output-rec\nExiting..." << endl;
      return 1;
    }

    //check files
    vector<string> checkfiles; //list of files to check for existence.
    if (vm.count("output-hap"))
    {
      cout << "Output corrected haps to:\t" << corrected_out << endl;
      checkfiles.push_back(corrected_out + ".haps");
      checkfiles.push_back(corrected_out + ".sample");
    }
    if (vm.count("output-generr"))
    {
      cout << "Output genotyping errors to:\t" << generr_out << endl;
      checkfiles.push_back(generr_out);
    }
    if (vm.count("output-rec"))
    {
      cout << "Output recombination map to:\t" << rec_out << endl;
      checkfiles.push_back(rec_out);
    }

    for (size_t i = 0; i < checkfiles.size(); i++)
    {
      if (fileexists(checkfiles[i]))
      {
        cerr << "ERROR: " << checkfiles[i] << " exists.  Will not overwrite.\nExiting..." << endl;
        return 1;
      }
    }
    if (gm.compare("") == 0)
      cout << "WARNING: You have not specified a genetic map! This is not recommended.\nRecombination rate will be set to 1.19 cM/Mb. (sex averaged)" << endl;

    //    pedhap ph(haps,fam,gm);
    pedhap ph(haps, gm, niteration);
    cout << "Correcting haplotypes based on pedigree structure..." << endl;
    ph.correct();

    if (vm.count("output-hap"))
      ph.haps->writeHaps(corrected_out);
    if (vm.count("output-generr"))
      ph.genotypingError(generr_out);
    if (vm.count("output-rec"))
      ph.recombinationMap(rec_out);
  }
  catch (exception &e)
  {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch (...)
  {
    cerr << "Exception of unknown type!\n";
  }

  return 0;
}
