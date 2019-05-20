#include "makePlots.h"
#include "PlotSetting.h"
#include <fstream>
#include <iostream>

int main(int argc, char** argv){
  
  int arg_tmp = -1;
  bool Is_TB = false;
  if(argc == 2){
    std::string tmp_str(argv[1]);
    arg_tmp = (int)argv[1][0] - int('0');
    
    if(arg_tmp == 0){
      Is_TB = true;
      cout << "Run TB Module\n" << endl;}
    if(arg_tmp == 1){
      cout << "Run NTU Bare PCB\n" << endl;}  }
  else if(argc == 1){
    Is_TB = true;
    cout << "Run TB Module\n" << endl;}
  else{
    cout << "wrong argument! EXIT!!" << endl;
    exit(-1);
  }

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
  M.Is_TB = Is_TB; 
  M.input_RUN = filename;
  M.Init();
  //M.PlotProducer();
  //M.Evt_display();
  M.Inj_Pulse_display();
  //  M.IdentifyInjCh();

		      
  
  return(0);
}
