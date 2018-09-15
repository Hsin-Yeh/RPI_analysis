#include "TApplication.h"
#include "TROOT.h"
#include "TChain.h"
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "makePlots.h"


int main (int argc, char* argv[])
{
  TChain *chain = new TChain("HGCBD");
  string filename;
  ifstream infile("input.txt");
  infile >> filename;
  infile.close();
  if( filename.length() > 2){
    chain->Add(filename.c_str());}
  else
    cout << "There is no input root file written in the input.txt!" << endl;

  std::cout << "To run it you type: " << std::endl;
  std::cout << "./makePlots"<< std::endl;

  makePlots makePlots(chain);
  makePlots.input_RUN=filename;
  makePlots.Loop();

  cout << "Last line in main ..." << endl;
  return 0;
}

