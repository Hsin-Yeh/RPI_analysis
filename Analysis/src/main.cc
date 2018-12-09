#include "makePlots.h"
#include "PlotSetting.h"
#include <fstream>
#include <iostream>

int main(){
  TChain *chain = new TChain("treeproducer/sk2cms");
  string filename;
  ifstream infile("data_input.txt");
  infile >> filename;
  infile.close();
  if( filename.length() > 2){
    chain->Add(filename.c_str());}
  else
    cout << "There is no input root file written in the input.txt!" << endl;

  makePlots M(chain);
  M.input_RUN = filename;
  M.Init();
  M.GainFactorReader();
  //M.Gain_factor_producer();
  //M.PlotProducer();
  //M.Evt_display();
  // M.Inj_Pulse_display();
  //  M.IdentifyInjCh();

		      
  
  return(0);
}
