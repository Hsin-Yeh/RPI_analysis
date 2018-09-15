#include "makePlots.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TImage.h"
#include <string>

ClassImp(hit)
ClassImp(hitcollection)
//Constructo
makePlots::makePlots(TChain* inchain):fChain(inchain)
{
  //  HITS = new hitcollection;
  readmap();
  cout << "Constructor of makePlot ... \n\n" << endl;
}
//Destructor
 makePlots::~makePlots()
{
  delete HITS;
  cout << "\n\n";
  cout << "Destructor of makePlot ... " << endl;
}

void makePlots::Init(){
  fChain->SetBranchAddress("hits",&HITS);
  outfile = new TFile("output.root","RECREATE");
}

void makePlots::Loop(){

  Init();
  app = new TApplication("app",0,0);
  TCanvas *c1 = new TCanvas;

  P_and_N(0,1);
  
  P_and_N(1,1);

}
  
Int_t makePlots::Cut(Long64_t entry, Long64_t sigma)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  if(entry<10*sigma){
    return -1;
  }
  else return 1;
}

void makePlots::P_and_N(int option,bool output){

  int skip_TS;
  int mem_of_SCA[NSCA][NCH][NSCA];
  int mem_of_ch [NSCA][NCH];
  int mem_of_chip[NSCA];
  int mem_of_module=0;
  int sumHG=0;
  int sumLG=0;
  double sumsqHG=0;
  double sumsqLG=0;
  int avg_HG=0;
  int avg_LG=0;
  int sigma_HG=0;
  int sigma_LG=0;
  int sumHG_SCA     [NCHIP][NCH][NSCA];
  int sumLG_SCA     [NCHIP][NCH][NSCA];
  double sumsqHG_SCA   [NCHIP][NCH][NSCA];
  double sumsqLG_SCA   [NCHIP][NCH][NSCA];
  int sumHG_ch  [NCHIP][NCH];
  int sumLG_ch  [NCHIP][NCH];
  double sumsqHG_ch [NCHIP][NCH];
  double sumsqLG_ch [NCHIP][NCH];
  int sumHG_chip [NCHIP];
  int sumLG_chip [NCHIP];
  double sumsqHG_chip [NCHIP];
  double sumsqLG_chip [NCHIP];
  double avg_HG_ch  [NCHIP][NCH];
  double sigma_HG_ch[NCHIP][NCH];
  double avg_LG_ch  [NCHIP][NCH];
  double sigma_LG_ch[NCHIP][NCH];
  double avg_HG_chip  [NCHIP];
  double sigma_HG_chip[NCHIP];
  double avg_LG_chip  [NCHIP];
  double sigma_LG_chip[NCHIP];

  double noisy_SCA_check[NCHIP][NCH][NSCA];
  double noisy_ch_check[NCHIP][NCH];

  char  name[50], title[50], filepath[50];
  
  TH1D* h_HGped[NCHIP][NCH][NSCA];
  TH1D* h_LGped[NCHIP][NCH][NSCA];
  TH1D* h_HGped_ch[NCHIP][NCH];
  TH1D* h_LGped_ch[NCHIP][NCH];
  TH1D* h_HGped_chip[NCHIP];
  TH1D* h_LGped_chip[NCHIP];
  TH1D* h_noisy_SCA = new TH1D("h_noisy_SCA","Noisy SCA Counts",10,0,10);
  TH1D* h_noisy_ch = new TH1D("h_noisy_ch","Noisy SCA Counts",10,0,10);

  TCanvas *c1 = new TCanvas();
  

  ////////// Initialize //////////
  for(int i = 0; i < NCHIP; ++i){
    sumHG_chip [i] = 0;
    sumLG_chip [i] = 0;
    sumsqHG_chip[i] = 0;
    sumsqLG_chip[i] = 0;
    mem_of_chip[i] = 0;
    
    sprintf(name,"h_HGpedChip%d",i);
    sprintf(title,"Pedestal_HG_Chip%d",i);
    h_HGped_chip[i] = new TH1D(name,title,500,-500,500);
    sprintf(name,"h_LGpedChip%d",i);
    sprintf(title,"Pedestal_LG_Chip%d",i);
    h_LGped_chip[i] = new TH1D(name,title,100,-200,200);

    for(int j = 0; j < NCH; ++j){
      
      sumHG_ch [i][j] = 0;
      sumLG_ch [i][j] = 0;
      sumsqHG_ch[i][j] = 0;
      sumsqLG_ch[i][j] = 0;
      mem_of_ch[i][j] = 0;

      sprintf(name,"h_HGpedChip%d_Ch%d",i,j);
      sprintf(title,"Pedestal_HG_Chip%d,Ch%d",i,j);
      h_HGped_ch[i][j] = new TH1D(name,title,100,-200,200);
      sprintf(name,"h_LGpedChip%d_Ch%d",i,j);
      sprintf(title,"Pedestal_LG_Chip%d,Ch%d",i,j);	
      h_LGped_ch[i][j] = new TH1D(name,title,100,-200,200);
      
      for(int k = 0; k < NSCA; ++k){
	
	sumHG_SCA     [i][j][k] = 0;
	sumsqHG_SCA   [i][j][k] = 0;
	sumLG_SCA     [i][j][k] = 0;
	sumsqLG_SCA   [i][j][k] = 0;
	mem_of_SCA[i][j][k] = 0;
	
	sprintf(name,"h_HGpedChip%d_Ch%d_SCA%d",i,j,k);
	sprintf(title,"Pedestal_HG_Chip%d,Ch%d_SCA%d",i,j,k);
	h_HGped[i][j][k] = new TH1D(name,title,100,-200,200);
	sprintf(name,"h_LGpedChip%d_Ch%d_SCA%d",i,j,k);
	sprintf(title,"Pedestal_LG_Chip%d,Ch%d_SCA%d",i,j,k);	
	h_LGped[i][j][k] = new TH1D(name,title,100,-200,200);
      }
    }
  }

  if(option == 0){
    cout << "calculate pedestal and noise based on TS 0 ~ 8" << endl;
    skip_TS = 9;  }
  
  else if(option == 1){
    cout << "calculate pedestal and noise based on TS 0 and 1"
	 << "(The test beam method)" << endl;
    skip_TS = 2;  }

  else{
    cout << "invalid option for pedestal noise calculation!" << endl;
    return;}

  int nevents = fChain->GetEntries();
  for(int ev = 0; ev < nevents ; ++ev){
    fChain -> GetEntry(ev);
    for(int i = 0 ; i < NSCA ; ++i)
      TS[i] = HITS->rollposition[i];

    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){
      H = HITS->Hits.at(hit); // H = hit (every channel would have one hit in an event)
      for(int sca = 0; sca < NSCA; ++sca){
	if(TS[sca] >= skip_TS) continue;
	sumHG += H.SCA_hg[sca];
	sumLG += H.SCA_lg[sca];
	sumsqHG += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumsqLG += H.SCA_lg[sca]*H.SCA_lg[sca];
	sumHG_chip [H.chip] += H.SCA_hg[sca];
	sumsqHG_chip[H.chip] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG_chip  [H.chip] += H.SCA_lg[sca];
	sumsqLG_chip[H.chip] += H.SCA_lg[sca]*H.SCA_lg[sca];
	sumHG_ch  [H.chip][H.ch] += H.SCA_hg[sca];
	sumsqHG_ch[H.chip][H.ch] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG_ch  [H.chip][H.ch] += H.SCA_lg[sca];
	sumsqLG_ch[H.chip][H.ch] += H.SCA_lg[sca]*H.SCA_lg[sca];
	sumHG_SCA  [H.chip][H.ch][sca] += H.SCA_hg[sca];
	sumsqHG_SCA[H.chip][H.ch][sca] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG_SCA  [H.chip][H.ch][sca] += H.SCA_lg[sca];
	sumsqLG_SCA[H.chip][H.ch][sca] += H.SCA_lg[sca]*H.SCA_lg[sca];
	mem_of_SCA[H.chip][H.ch][sca]++;
	mem_of_ch[H.chip][H.ch]++;
	mem_of_chip[H.chip]++;
	mem_of_module++;
      }
    }
  }
  
  avg_HG = (double)sumHG/mem_of_module;
  sigma_HG=sigmaCal(mem_of_module,sumHG,sumsqHG);
  avg_LG = (double)sumLG/mem_of_module;
  sigma_LG=sigmaCal(mem_of_module,sumLG,sumsqLG);
  for(int i = 0; i < NCHIP; ++i){
    avg_HG_chip  [i] = (double)sumHG_chip[i]/mem_of_chip[i];
    sigma_HG_chip[i]=sigmaCal(mem_of_chip[i],sumHG_chip[i],sumsqHG_chip[i]);
    avg_LG_chip  [i] = (double)sumLG_chip[i]/mem_of_chip[i];
    sigma_LG_chip[i]=sigmaCal(mem_of_chip[i],sumLG_chip[i],sumsqLG_chip[i]);
    for(int j = 0; j < NCH; ++j){
      avg_HG_ch  [i][j] = (float)sumHG_ch[i][j]/mem_of_ch[i][j];
      sigma_HG_ch[i][j]=sigmaCal(mem_of_ch[i][j],sumHG_ch[i][j],sumsqHG_ch[i][j]);
      avg_LG_ch  [i][j] = (float)sumLG_ch[i][j]/mem_of_ch[i][j];
      sigma_LG_ch[i][j]=sigmaCal(mem_of_ch[i][j],sumLG_ch[i][j],sumsqLG_ch[i][j]);
      noisy_ch_check[i][j] = (sigma_HG_ch[i][j] - sigma_HG_chip[i])/sigma_HG_chip[i];
      h_noisy_ch->Fill(noisy_ch_check[i][j]);
      for(int k = 0; k < NSCA; ++k){
	avg_HG_SCA  [i][j][k] = (float)sumHG_SCA  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_HG_SCA[i][j][k]=sigmaCal(mem_of_SCA[i][j][k],sumHG_SCA[i][j][k],sumsqHG_SCA[i][j][k]);
	avg_LG_SCA  [i][j][k] = (float)sumLG_SCA  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_LG_SCA[i][j][k]=sigmaCal(mem_of_SCA[i][j][k],sumLG_SCA[i][j][k],sumsqLG_SCA[i][j][k]);
	noisy_SCA_check[i][j][k] = (sigma_HG_SCA[i][j][k] - sigma_HG_ch[i][j])/sigma_HG_ch[i][j];
	h_noisy_SCA->Fill(noisy_SCA_check[i][j][k]);
	if(noisy_SCA_check[i][j][k]>1){
	  //cout << "chip_" << i << "ch_" << j << "SCA" << k << endl;
	}
      }
      if(noisy_ch_check[i][j]>1){
	//cout << "chip_" << i << "ch_" << j << endl;
      }
      delete h_HGped_ch[i][j];
      delete h_LGped_ch[i][j];
    }
  }

  
  if(output == true){
    string filename;
    filename = input_RUN;
    int start = filename.find_last_of("/");
    int end   = filename.find(".root");
    string outf = filename.substr(start+1,end-start-1);
    char outtitleH[100];
    char outtitleL[100];   
    sprintf(outtitleH,"long_term_ped/opt%d/%s_HG.txt",option,outf.c_str());
    ofstream fileHG(outtitleH);
    sprintf(outtitleL,"long_term_ped/opt%d/%s_LG.txt",option,outf.c_str());
    ofstream fileLG(outtitleL);

    fileHG << fixed << setprecision(2) << avg_HG << " ";
    fileLG << fixed << setprecision(2) << avg_LG << " ";
    fileHG << "\n";
    fileLG << "\n";
    fileHG << fixed << setprecision(2) << sigma_HG << " ";
    fileLG << fixed << setprecision(2) << sigma_LG << " ";
    fileHG << "\n";
    fileLG << "\n";
    fileHG << "CHIP";
    fileLG << "CHIP";
    fileHG << "\n";
    fileLG << "\n";
    for(int i = 0; i < NCHIP; i++){
      fileHG << i << "\t" ;
      fileLG << i << "\t" ;
      fileHG << fixed << setprecision(2) << avg_HG_chip[i] << " ";
      fileLG << fixed << setprecision(2) << avg_LG_chip[i] << " ";
      fileHG << "\n";
      fileLG << "\n";
      fileHG << i << "\t" ;
      fileLG << i << "\t" ;
      fileHG << fixed << setprecision(2) << sigma_HG_chip[i] << " ";
      fileLG << fixed << setprecision(2) << sigma_LG_chip[i] << " ";
      fileHG << "\n";
      fileLG << "\n";
    }
    fileHG << "CHIP\tCH";
    fileLG << "CHIP\tCH";
    fileHG << "\n";
    fileLG << "\n";
    for(int i = 0; i < NCHIP; ++i){
      for(int j = 0; j < NCH; ++j){
	fileHG << i << "\t" << j << "\t";
	fileLG << i << "\t" << j << "\t";
	fileHG << fixed << setprecision(2) << avg_HG_ch[i][j] << " ";
	fileLG << fixed << setprecision(2) << avg_LG_ch[i][j] << " ";
	fileHG << "\n";
	fileLG << "\n";
	fileHG << i << "\t" << j << "\t";
	fileLG << i << "\t" << j << "\t";
	fileHG << fixed << setprecision(2) << sigma_HG_ch[i][j] << " ";
	fileLG << fixed << setprecision(2) << sigma_LG_ch[i][j] << " ";
	fileHG << "\n";
	fileLG << "\n";
      }
    }
    
    fileHG << "CHIP\tCH\t";
    fileLG << "CHIP\tCH\t";
    for(int i = 0; i < NSCA; ++i){
      fileHG << "SCA " << i << " ";
      fileLG << "SCA " << i << " ";}
    fileHG << "\n";
    fileLG << "\n";
    for(int i = 0; i < NCHIP; ++i){
      for(int j = 0; j < NCH; ++j){
	fileHG << i << "\t" << j << "\t";
	fileLG << i << "\t" << j << "\t";
	for(int k = 0; k < NSCA; ++k){
	  fileHG << fixed << setprecision(2) << avg_HG_SCA[i][j][k] << " ";
	  fileLG << fixed << setprecision(2) << avg_LG_SCA[i][j][k] << " ";}
	fileHG << "\n";
	fileLG << "\n";
	fileHG << i << "\t" << j << "\t";
	fileLG << i << "\t" << j << "\t";

	for(int k = 0; k < NSCA; ++k){
	  fileHG << fixed << setprecision(2) << sigma_HG_SCA[i][j][k] << " ";
	  fileLG << fixed << setprecision(2) << sigma_LG_SCA[i][j][k] << " ";}
	fileHG << "\n";
	fileLG << "\n";      
      }
    }
    fileHG.close();
    fileLG.close();
    cout << "output mode is selected output file will be:" << endl;
    cout << "1. " << outtitleH << "\n" << "2. " << outtitleL << endl;
  }

}


void makePlots::readmap(){
  ifstream file("./src_txtfile/CH_map.txt");
  string line;
  int chip,ch,type,formatCH;
  double posx,posy;
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> chip >> ch >> posx >> posy >> type;
    formatCH = chip*32+ch/2;
    CHmap[formatCH] = make_pair(posx,posy);}
  file.close();
  //Since there is no such pad, assign a unreasonable value
  CHmap[2*32+60/2] = make_pair(1000.,1000.);

}

double makePlots::sigmaCal(int N, int sum, double sqsum){
  double avgsum = sum/N;
  double avgsqsum = sqsum/N;
  avgsqsum -= avgsum*avgsum;
  avgsqsum = sqrt(avgsqsum);
  
  return avgsqsum;
  
}

void makePlots::InitTH2Poly(TH2Poly& poly)
{
  int MAXVERTICES = 6;
  double HexX[MAXVERTICES];
  double HexY[MAXVERTICES];
  int iu,iv,CellXYsize;
  ifstream file("src_txtfile/poly_frame.txt");
  string line;

  
  for(int header = 0; header < 4; ++header )     getline(file,line);
  
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> iu >> iv >> CellXYsize;    
    for(int i = 0; i < CellXYsize ; ++i){
      getline(file,line);
      file >> HexX[i] >> HexY[i];
    }
    
    /*if(CellXYsize == 4){
      CellXYsize++;
      HexX[4] = HexX[0];
      HexY[4] = HexY[0];}*/
    poly.AddBin(CellXYsize, HexX, HexY);
  }
  file.close();

}
