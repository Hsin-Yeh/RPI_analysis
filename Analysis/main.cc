#include "makePlots.h"
#include <fstream>
#include <iostream>

int main(){
  TChain *chain = new TChain("HGCBD");
  string filename;
  ifstream infile("input.txt");
  infile >> filename;
  infile.close();
  if( filename.length() > 2){
    chain->Add(filename.c_str());}
  else
    cout << "There is no input root file written in the input.txt!" << endl;

  makePlots M(chain);
  M.input_RUN = filename;
  M.Init();
  //M.Loop();
  M.Evt_display();
  //M.Inj_Pulse_display();
  
  return(0);
}
