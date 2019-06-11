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
int displayChannel = -1;
int acq_type = 0;
int lowerR = -1, upperR = -1;
bool subPed_flag = false;
bool maskCh_flag = false;
bool pulseDisplay_flag = false;


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
	  pulseDisplay_flag = true;
	  if ( isNumber( arg_list[iarg+1] ) ) {
		displayChannel = atoi(arg_list[iarg+1].c_str());
		iarg+=2;
	  }
	  else iarg++;
	}
	else if ( arg_list[iarg] == "-i" ) {
	  pulseDisplay_flag = true;
	  acq_type = 1; 
	  iarg++;
	}
	else if ( arg_list[iarg] == "-s" ) {
	  pulseDisplay_flag = true;
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
	else {
	  std::cout << "Unknown option... print usage" << std::endl;
	}
  }
	
  if ( pulseDisplay_flag )
	main_plotPulseDisplay();
  else
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
  M.PlotProducer();
}

void main_plotPulseDisplay() {
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
  M.Pulse_display( displayChannel, acq_type, lowerR, upperR );
		      
}

 bool isNumber(string s) { 
   for (int i = 0; i < s.length(); i++) 
	 if (isdigit(s[i]) == false) 
	   return false; 

   return true; 
 } 

