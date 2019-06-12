#include "makePlots.h"
#include "PlotSetting.h"
#include <fstream>
#include <iostream>


// default Setup
string main_outpath = "./";
string main_datainput = "./data_input.txt";
string pedfile  = "./pedestal";
string gainfile = "./src_txtfile/TPro_fittingoutput.txt";

// Utility
void main_plotPulseDisplay(); 
void main_makePlots();
bool isNumber(string s);
int  displayChannel = -1;
int  acq_type = 0;
int  lowerR = -1, upperR = -1;
bool subPed_flag = true;
bool maskCh_flag = false;
int  anaType = 0;


int main(int argc, char** argv){
  

  string arg_string;
  vector<string> arg_list;
  for(int i = 0 ; i < argc ; ++i){
	arg_string = argv[i];
	arg_list.push_back(arg_string);
  }

  int iarg = 1;
  while ( iarg < argc ) {
	if ( arg_list[iarg] == "-p" ) {
	  anaType = 1;
	  cout << arg_list[iarg+1] << endl;
	  if ( isNumber( arg_list[iarg+1] ) ) {
		displayChannel = atoi(arg_list[iarg+1].c_str());
		iarg+=2;
	  }
	  else iarg++;
	}
	else if ( arg_list[iarg] == "-i" ) {
	  anaType = 1;
	  acq_type = 1; 
	  iarg++;
	}
	else if ( arg_list[iarg] == "-s" ) {
	  anaType = 1;
	  acq_type = 2;
	  iarg++;
	}
	else if ( arg_list[iarg] == "-r" ) {
	  lowerR =  atoi(arg_list[iarg+1].c_str());
	  upperR =  atoi(arg_list[iarg+2].c_str());
	  iarg+=3;
	}
	else if ( arg_list[iarg] == "-n" || arg_list[iarg] == "-noSubPed" ) {
	  subPed_flag = false;
	  iarg++;
	}
	else if ( arg_list[iarg] == "-m" || arg_list[iarg] == "-mask" ) {
	  maskCh_flag = true;
	  iarg++;
	}
	else if ( arg_list[iarg] == "-c" || arg_list[iarg] == "-cosmic" ) {
	  anaType = 2;
	  iarg++;
	}
	else {
	  std::cout << "Unknown option... print usage" << std::endl;
	}
  }
  main_makePlots();
  
  return (0);
}

void main_makePlots() {
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
  M.input_fileName = filename;
  M.maskCh_flag = maskCh_flag;
  M.subPed_flag = subPed_flag;
  M.Init( pedfile, gainfile );
  if ( anaType == 0 )
	M.PlotProducer();
  else if ( anaType == 1 )
	M.Pulse_display();
  else
	M.cosmicAnalyzer();
}

bool isNumber(string s) { 
   for (int i = 0; i < s.length(); i++) 
	 if (isdigit(s[i]) == true) 
	   return true; 

   return false; 
 } 

