#include "makePlots.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TImage.h"

ClassImp(hit)
ClassImp(hitcollection)
//Constructo
makePlots::makePlots(TChain* inchain):fChain(inchain)
{
  //HITS = new hitcollection;
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
  app = new TApplication("app",0,0);
  c1 = new TCanvas();
}

void makePlots::Loop(){

  //P_and_N(0,1);
  //read_P_and_N("ped_result/Module_1_RUN_300318_0527");


  //==================== Call the Parameters ====================

  char title[50];
  char plot_title[50];
  int nevents = fChain->GetEntries();
  int injch;
  int injADC=0;
  int injevents_perdac = 1;
  int injevents = nevents/injevents_perdac;

  TH1D *h = new TH1D("h","",100,150,250); //("title","",slice,star,end)
  TH1D *h_TOTS = new TH1D("h_TOTS","",100,5,500);
  TH1D *h_TOTF = new TH1D("h_TOTF","",100,1000,3000);
  TH1D *h_TOAR = new TH1D("h_TOAR","",100,1000,3000);
  TH1D *h_TOAF = new TH1D("h_TOAF","",100,1000,3000);
  
  int ADC_H[injevents], ADC_L[injevents], TOTS[injevents], dac_ctrl[injevents];
  int ADC_event[nevents], ADC[injevents_perdac], n[injevents_perdac];
  int Crosstalk_ADC_H[6][injevents], Crosstalk_ADC_L[6][injevents], Crosstalk_TOTS[6][injevents];
  int NoisyChannel_ADC_H[injevents];;
  
  int dac = -1; int test =0;
  int cross_num = 6;
  char leg[50], img_title[50];
  TMultiGraph* mg = new TMultiGraph();
  TGraph **g = new TGraph*[13];
  TLegend *legend = new TLegend(0.85,0.8,1.,1.);
  //legend->SetNColumns(2);
  TImage *img = TImage::Create();

  // Set Output Root File
  string Input_filename;
  Input_filename = input_RUN;
  int start = Input_filename.find_last_of("/");
  int end   = Input_filename.find(".root");
  string outf = Input_filename.substr(start+1,end-start-1);
  sprintf(title,"output_root/%s.root",outf.c_str());
  //  TFile outf(title,"recreate");
  TTree *outT1 = new TTree("Rechit_var","Rechit_var");
  
  
  //==================== Initialize ====================
  
  //for(int ts = 0; ts < 7; ts++){
  int ts = 4; //choose this time sample to be the peak
  dac = -1;

  for (int i=0; i<injevents; i++){
    ADC_H[i] = 0;
    ADC_L[i] = 0;
    ADC_event[i] = 0;
    dac_ctrl[i] = i;
  }

    //==================== Loop over the events ====================
   
  for(int ev = 0; ev < nevents ; ++ev){  
    fChain -> GetEntry(ev); //== Get HITcollection from root file event 

    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITS->rollposition[i];
    
    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){
      H = HITS->Hits.at(hit);
      if(!H.CCorNC) continue; // remove unconnected channel
      if(H.ch!=HITS->inj_ch.at(0)) continue;
      for(int sca=0; sca<NSCA; sca++){
	if(TS[sca]==3){
	  ADC_H[ev] = H.SCA_hg[sca];
	  ADC_L[ev] = H.SCA_lg[sca];
	}
      }
    }
  }
  
    //==================== End of Loop ====================
    //
    //==================== Draw Plots ====================
    
    //g[ts] = new TGraph(injevents,dac_ctrl,ADC_H);
    TGraph* gh = new TGraph(nevents,dac_ctrl,ADC_H);
    TGraph* gl = new TGraph(nevents,dac_ctrl,ADC_L);
    TGraph* gTOT = new TGraph(injevents,dac_ctrl,TOTS);
    TGraph* ghlratio = new TGraph(injevents,ADC_L,ADC_H);
    TGraph* gltratio = new TGraph(injevents,TOTS,ADC_L);
    TGraph** gcross_h = new TGraph*[cross_num];
    TGraph** gcross_l = new TGraph*[cross_num];
    TGraph** gcross_TOTS = new TGraph*[cross_num];
    TGraph* gnoisy_h = new TGraph(injevents,dac_ctrl,NoisyChannel_ADC_H);
    
    
    gh->GetYaxis()->SetRangeUser(0,200);
    //    gl->GetXaxis()->SetRangeUser(xmin,xmax);
    sprintf(plot_title,"High Gain");
    gh->SetTitle(plot_title);
    gh->GetXaxis()->SetTitle("Event");
    gh->GetYaxis()->SetTitle("ADC");
    gh->SetMarkerStyle(7);
    gh->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();
    
    sprintf(plot_title,"Low Gain");
    gl->SetTitle(plot_title);
    gl->GetXaxis()->SetTitle("Event");
    gl->GetYaxis()->SetTitle("ADC");
    gl->SetMarkerStyle(7);
    gl->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();

  //=================== End of filling hist =======================
  
}

void makePlots::Evt_display(){

  int nevents = fChain->GetEntries();
  char plot_title[50];
  for(int ev = 0; ev < nevents ; ++ev){
    if(ev % 1 != 0) continue;

    fChain -> GetEntry(ev); // Get HITcollection from root file event
    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);
    
    poly->Draw("colztext0");
  
    sprintf(plot_title,"test%d",ev);
    poly->SetTitle(plot_title);
    c1->Update();
    gPad->WaitPrimitive();
        
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITS->rollposition[i];

    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){  
      H = HITS->Hits.at(hit);
      if(!H.CCorNC) continue; // remove unconnected channel
      int forCH = H.chip*32+H.ch/2;
      float X = CHmap[forCH].first;
      float Y = CHmap[forCH].second;	  

      for(int sca = 0; sca < NSCA; ++sca){
	if(TS[sca] == 8 ){
	  cout << H.SCA_lg[sca] << endl;
	  //	  H.SCA_hg[sca] -= avg_HG[H.chip][H.ch][sca]; // pedestal subtraction
	  //poly->Fill(X,Y,H.SCA_lg[sca]);
	  poly->Fill(X,Y,H.ch);
	}
      }
    }  
    
    poly->Draw("colztext0");
  
    sprintf(plot_title,"HG_TS4_evt%d",ev);
    poly->SetTitle(plot_title);
    c1->Update();
    gPad->WaitPrimitive();
    delete poly;
  }
}



void makePlots::Inj_Pulse_display(){

  TGraph *gr;
  int nevents = fChain->GetEntries();
  char plot_title[50];
  
  for(int ev = 0; ev < nevents ; ++ev){
    if(ev % 35 != 0) continue;
    fChain -> GetEntry(ev);
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITS->rollposition[i];    
    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){
      H = HITS->Hits.at(hit);
      if(H.ch != HITS->inj_ch.at(0)) continue;
      for(int sca=0; sca<13; sca++){
	if(TS[sca]==8) cout << H.SCA_lg[sca] << endl;
      }

      gr = new TGraph(13, TS,H.SCA_lg );
      gr->SetMarkerColor(H.chip+2);
      gr->SetMarkerStyle(22);
      gr->SetMarkerSize(1.2);
      gr->Draw("AP");
      sprintf(plot_title,"HG_evt %d chip %d",ev,H.chip);
      gr->SetTitle(plot_title);
      gr->GetXaxis()->SetTitle("TS");
      gr->GetYaxis()->SetTitle("ADC");

      c1->Update();
      //if(ev == 400){
      //sprintf(plot_title,"%s.png",plot_title);
      //c1->SaveAs(plot_title);} // remove the comment to save plots
      gPad->WaitPrimitive();

    }
  }
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




void makePlots::Crosstalk(Int_t CH){
  int chip = 3;
  int cross_num = 6;
  int forCH = chip*32+CH/2;
  float Xdist = 0.974452;
  float Ydist = 0.5626;
  float cross_posx[cross_num], cross_posy[cross_num];
  float X = CHmap[forCH].first;
  float Y = CHmap[forCH].second;
  
  cross_posx[0] = X - Xdist;
  cross_posy[0] = Y - Ydist;
  cross_posx[1] = X - Xdist;
  cross_posy[1] = Y + Ydist;
  cross_posx[2] = X;
  cross_posy[2] = Y - 2*Ydist;
  cross_posx[3] = X;
  cross_posy[3] = Y + 2*Ydist;
  cross_posx[4] = X + Xdist;
  cross_posy[4] = Y - Ydist;
  cross_posx[5] = X + Xdist;
  cross_posy[5] = Y + Ydist;

  for(int i=0; i<6; i++){
    //cout << cross_posx[i] << " " << cross_posy[i] << endl;
    int ch = 0;
    
    while(abs(CHmap[ch].first-cross_posx[i]) > 1e-4 || abs(CHmap[ch].second-cross_posy[i]) > 1e-4){
      cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
      ch++;
      if(ch>256) break;
    }
    //    cross_ch[i] = (ch-chip*32)*2;
    cross_ch[i] = ch;
  }
  for(int i=0; i<6; i++){
    //cout << cross_ch[i] << endl;
  }
  
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
	  cout << "chip_" << i << "ch_" << j << "SCA" << k << endl;
	}
      }
      if(noisy_ch_check[i][j]>1){
	cout << "chip_" << i << "ch_" << j << endl;
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
    sprintf(outtitleH,"long_term_ped/%s_HG.txt",outf.c_str());
    ofstream fileHG(outtitleH);
    sprintf(outtitleL,"long_term_ped/%s_LG.txt",outf.c_str());
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

void makePlots::Injection_ana(int Inj_ch){

  TCanvas *c1 = new TCanvas();

  char plot_title[50];
  int nevents = fChain->GetEntries();
  int injevents_perdac = 1;
  int injevents = nevents/injevents_perdac;

  TH1D *h = new TH1D("h","",100,150,250); //("title","",slice,star,end)
  TH1D *h_TOTS = new TH1D("h_TOTS","",100,5,500);
  TH1D *h_TOTF = new TH1D("h_TOTF","",100,1000,3000);
  TH1D *h_TOAR = new TH1D("h_TOAR","",100,1000,3000);
  TH1D *h_TOAF = new TH1D("h_TOAF","",100,1000,3000);
  
  int Inj_hg[injevents], Inj_lg[injevents], Inj_tots[injevents], Inj_toar[injevents], Inj_toaf[injevents];
  int ADC_event[nevents], ADC[injevents_perdac], n[injevents_perdac], dac_ctrl[injevents];
  int Crosstalk_ADC_H[6][injevents], Crosstalk_ADC_L[6][injevents], Crosstalk_TOTS[6][injevents];
  int NoisyChannel_ADC_H[injevents];;
  
  int dac = -1; int test =0;
  int cross_num = 6;
  int ts = 4;
  char leg[50], img_title[50];
  TMultiGraph* mg = new TMultiGraph();
  TGraph **g = new TGraph*[13];
  TLegend *legend = new TLegend(0.85,0.8,1.,1.);
  //  legend->SetNColumns(2);
  TImage *img = TImage::Create();

  for(int i=0; i < injevents; i++){
    dac_ctrl[i] = i;
  }
  
  //==================== Loop over the events ==================== 

  
  for(int ev = 0; ev < nevents ; ++ev){  
    fChain -> GetEntry(ev); //== Get HITcollection from root file event 

    // Update the dac number
    if(ev%injevents_perdac == 0){
      dac++; 
      //sprintf(plot_title,"DAC:%d",dac_ctrl[dac]);
      //TGraph *g = new TGraph(injevents_perdac,n,ADC);
      //g->SetTitle(plot_title);
      //g->Draw("A*");
      //c1->Update();
      //gPad->WaitPrimitive();
      delete g;
    }
       
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITS->rollposition[i];
    double HGSubPed, LGSubPed;
    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){
      H = HITS->Hits.at(hit);
      for(int sca = 0; sca < NSCA; ++sca){
	if( TS[sca] == ts ){	    
	  HGSubPed = H.SCA_hg[sca] - avg_HG_SCA[H.chip][H.ch][sca]; // pedestal subtraction
	  LGSubPed = H.SCA_lg[sca] - avg_LG_SCA[H.chip][H.ch][sca];
	  if (H.ch == Inj_ch){
	    Inj_hg[dac] = HGSubPed;
	    Inj_lg[dac] = LGSubPed;
	    Inj_tots[dac] = H.TOTS;
	    Inj_toar[dac] = H.TOAR;
	    Inj_toaf[dac] = H.TOAF;
	  }
	  else {
	  }
	}
      }
    }
  }
  
  
    for(int i=0; i<injevents; i++){
      Inj_hg[i] /= injevents_perdac;
      Inj_lg[i] /= injevents_perdac;
    }

    TGraph* g_Inj_hg = new TGraph(injevents,dac_ctrl,Inj_hg);
    TGraph* g_Inj_lg = new TGraph(injevents,dac_ctrl,Inj_lg);
    TGraph* g_Inj_tots = new TGraph(injevents,dac_ctrl,Inj_tots);
    TGraph* g_hlratio = new TGraph(injevents,Inj_lg,Inj_hg);
    TGraph* g_ltratio = new TGraph(injevents,Inj_tots,Inj_lg);

    sprintf(plot_title,"High Gain");
    g_Inj_hg->SetTitle(plot_title);
    g_Inj_hg->GetXaxis()->SetTitle("Event");
    g_Inj_hg->GetYaxis()->SetTitle("ADC");
    g_Inj_hg->SetMarkerStyle(7);
    g_Inj_hg->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();
    c1->Write();
    
    sprintf(plot_title,"Low Gain");
    g_Inj_lg->SetTitle(plot_title);
    g_Inj_lg->GetXaxis()->SetTitle("Event");
    g_Inj_lg->GetYaxis()->SetTitle("ADC");
    g_Inj_lg->SetMarkerStyle(7);
    g_Inj_lg->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();
    c1->Write();

    sprintf(plot_title,"TOTS");
    g_Inj_tots->SetTitle(plot_title);
    g_Inj_tots->GetXaxis()->SetTitle("Event");
    g_Inj_tots->GetYaxis()->SetTitle("ADC");
    g_Inj_tots->SetMarkerStyle(7);
    g_Inj_tots->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();
    c1->Write();


    sprintf(plot_title,"High vs Low Gain");
    g_hlratio->SetTitle(plot_title);
    g_hlratio->GetXaxis()->SetTitle("Low");
    g_hlratio->GetYaxis()->SetTitle("High");
    g_hlratio->SetMarkerStyle(24);
    g_hlratio->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();
    c1->Write();

    sprintf(plot_title,"Low Gain vs TOT");
    g_ltratio->SetTitle(plot_title);
    g_ltratio->GetXaxis()->SetTitle("TOT");
    g_ltratio->GetYaxis()->SetTitle("Low");
    g_ltratio->SetMarkerStyle(24);
    g_ltratio->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();
    c1->Write();
    outfile->Write();

    //==================== End of Loop ====================

    delete c1;

}

void makePlots::Pedestal_ana(int option){
  int skip_TS;
  TFile* outfile2;
  outfile2 = new TFile("output2.root","RECREATE");

  
  int mem_of_SCA[NSCA][NCH][NSCA];
  int mem_of_ch [NSCA][NCH];
  int mem_of_chip[NSCA];
  int sumHG     [NCHIP][NCH][NSCA];
  int sumLG     [NCHIP][NCH][NSCA];
  double sumsqHG   [NCHIP][NCH][NSCA];
  double sumsqLG   [NCHIP][NCH][NSCA];
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
	
	sumHG     [i][j][k] = 0;
	sumsqHG   [i][j][k] = 0;
	sumLG     [i][j][k] = 0;
	sumsqLG   [i][j][k] = 0;
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
  

  ////////// Select method to calculate the pedestal //////////
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
  
  
  ////////// Loop over events //////////
  
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
	sumHG_chip [H.chip] += H.SCA_hg[sca];
	sumsqHG_chip[H.chip] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG_chip  [H.chip] += H.SCA_lg[sca];
	sumsqLG_chip[H.chip] += H.SCA_lg[sca]*H.SCA_lg[sca];
	sumHG_ch  [H.chip][H.ch] += H.SCA_hg[sca];
	sumsqHG_ch[H.chip][H.ch] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG_ch  [H.chip][H.ch] += H.SCA_lg[sca];
	sumsqLG_ch[H.chip][H.ch] += H.SCA_lg[sca]*H.SCA_lg[sca];
	sumHG  [H.chip][H.ch][sca] += H.SCA_hg[sca];
	sumsqHG[H.chip][H.ch][sca] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG  [H.chip][H.ch][sca] += H.SCA_lg[sca];
	sumsqLG[H.chip][H.ch][sca] += H.SCA_lg[sca]*H.SCA_lg[sca];
	mem_of_SCA[H.chip][H.ch][sca]++;
	mem_of_ch[H.chip][H.ch]++;
	mem_of_chip[H.chip]++;
	/*
	h_HGped[H.chip][H.ch][sca]->Fill(H.SCA_hg[sca]);
	h_LGped[H.chip][H.ch][sca]->Fill(H.SCA_lg[sca]);
	h_HGped_ch[H.chip][H.ch]->Fill(H.SCA_hg[sca]);
	h_LGped_ch[H.chip][H.ch]->Fill(H.SCA_lg[sca]);
	h_HGped_chip[H.chip]->Fill(H.SCA_hg[sca]);
	h_LGped_chip[H.chip]->Fill(H.SCA_lg[sca]);	
	*/
      }
    }
  }

  //////////  End of loop  //////////

  //////////  Draw Plots  //////////
  for(int i = 0; i < NCHIP; ++i){
    avg_HG_chip  [i] = (double)sumHG_chip[i]/mem_of_chip[i];
    avg_LG_chip  [i] = (double)sumLG_chip[i]/mem_of_chip[i];
    sumHG_chip[i] = 0;
    sumLG_chip[i] = 0;
    sumsqHG_chip[i] = 0;
    sumsqLG_chip[i] = 0;
    mem_of_chip[i] = 0;
    for(int j = 0; j < NCH; ++j){
      avg_HG_ch  [i][j] = (double)sumHG_ch[i][j]/mem_of_ch[i][j];
      avg_LG_ch  [i][j] = (double)sumLG_ch[i][j]/mem_of_ch[i][j];
      sumHG_ch[i][j] = 0;
      sumLG_ch[i][j] = 0;
      sumsqHG_ch[i][j] = 0;
      sumsqLG_ch[i][j] = 0;
      mem_of_ch[i][j] = 0;


      for(int k = 0; k < NSCA; ++k){
	avg_HG_SCA  [i][j][k] = (double)sumHG  [i][j][k]/mem_of_SCA[i][j][k];
	avg_LG_SCA  [i][j][k] = (double)sumLG  [i][j][k]/mem_of_SCA[i][j][k];
	sumHG[i][j][k] = 0;
	sumLG[i][j][k] = 0;
	sumsqHG[i][j][k] = 0;
	sumsqLG[i][j][k] = 0;

	mem_of_SCA[i][j][k] = 0;


      }
    }
  }

  ////////// Loop over events //////////
  
  for(int ev = 0; ev < nevents ; ++ev){
    fChain -> GetEntry(ev);
    for(int i = 0 ; i < NSCA ; ++i)
      TS[i] = HITS->rollposition[i];

    int nhits = HITS->hit_num;
    double HGSubPed, LGSubPed;
    for(int hit = 0; hit < nhits ; ++hit){
      H = HITS->Hits.at(hit); // H = hit (every channel would have one hit in an event)
      for(int sca = 0; sca < NSCA; ++sca){
	if(TS[sca] >= skip_TS) continue;
	HGSubPed = (double)H.SCA_hg[sca] - avg_HG_SCA[H.chip][H.ch][sca];
	LGSubPed = (double)H.SCA_lg[sca] - avg_LG_SCA[H.chip][H.ch][sca];
	sumHG_chip [H.chip] += HGSubPed;
	sumsqHG_chip[H.chip] += HGSubPed*HGSubPed;
	sumLG_chip  [H.chip] += LGSubPed;
	sumsqLG_chip[H.chip] += LGSubPed*LGSubPed;
	sumHG_ch  [H.chip][H.ch] += HGSubPed;
	sumsqHG_ch[H.chip][H.ch] += HGSubPed*HGSubPed;
	sumLG_ch  [H.chip][H.ch] += LGSubPed;
	sumsqLG_ch[H.chip][H.ch] += LGSubPed*LGSubPed;
	sumHG  [H.chip][H.ch][sca] += HGSubPed;
	sumsqHG[H.chip][H.ch][sca] += HGSubPed*HGSubPed;
	sumLG  [H.chip][H.ch][sca] += LGSubPed;
	sumsqLG[H.chip][H.ch][sca] += LGSubPed*LGSubPed;
	h_HGped[H.chip][H.ch][sca]->Fill(HGSubPed);
	h_LGped[H.chip][H.ch][sca]->Fill(LGSubPed);
	h_HGped_ch[H.chip][H.ch]->Fill(HGSubPed);
	h_LGped_ch[H.chip][H.ch]->Fill(LGSubPed);
	h_HGped_chip[H.chip]->Fill(HGSubPed);
	h_LGped_chip[H.chip]->Fill(LGSubPed);
	mem_of_SCA[H.chip][H.ch][sca]++;
	mem_of_ch[H.chip][H.ch]++;
	mem_of_chip[H.chip]++;
      }
    }
  }

  //////////  End of loop  //////////

  

  for(int i = 0; i < NCHIP; ++i){
    cout << sumHG_chip[i] << endl;
    avg_HG_chip  [i] = (double)sumHG_chip[i]/mem_of_chip[i];
    sigma_HG_chip[i]=sigmaCal(mem_of_chip[i],sumHG_chip[i],sumsqHG_chip[i]);
    avg_LG_chip  [i] = (double)sumLG_chip[i]/mem_of_chip[i];
    sigma_LG_chip[i]=sigmaCal(mem_of_chip[i],sumLG_chip[i],sumsqLG_chip[i]);
    cout << sigma_HG_chip[i] << endl;
    cout << mem_of_chip[i] << endl;
    cout << avg_HG_chip[i] << endl;
    //h_HGped_chip[i]->Draw();
    //c1->Update();
    //gPad->WaitPrimitive();
    sprintf(name,"h_HGped%d",i);
    sprintf(filepath,"../plots/Pedestal/Pedestal_HG_Chip%d.pdf",i);
    //c1->SetTitle(name);
    //c1->SaveAs(filepath);
    cout << sigma_LG_chip[i] << endl;
    cout << mem_of_chip[i] << endl;
    cout << avg_LG_chip[i] << endl;
    //h_LGped_chip[i]->Draw();
    //c1->Update();
    //gPad->WaitPrimitive();
    sprintf(name,"h_LGped%d",i);
    sprintf(filepath,"../plots/Pedestal/Pedestal_LG_Chip%d.pdf",i);
    //c1->SetTitle(name);
    //c1->SaveAs(filepath);

    
    for(int j = 0; j < NCH; ++j){
      avg_HG_ch  [i][j] = (float)sumHG_ch[i][j]/mem_of_ch[i][j];
      sigma_HG_ch[i][j]=sigmaCal(mem_of_ch[i][j],sumHG_ch[i][j],sumsqHG_ch[i][j]);
      avg_LG_ch  [i][j] = (float)sumLG_ch[i][j]/mem_of_ch[i][j];
      sigma_LG_ch[i][j]=sigmaCal(mem_of_ch[i][j],sumLG_ch[i][j],sumsqLG_ch[i][j]);
      noisy_ch_check[i][j] = (sigma_HG_ch[i][j] - sigma_HG_chip[i])/sigma_HG_chip[i];
      h_noisy_ch->Fill(noisy_ch_check[i][j]);
      cout << sigma_HG_ch[i][j] << endl;
      cout << avg_HG_ch[i][j] <<endl;
      cout << mem_of_ch[i][j] << endl;
      /*
      h_HGped_ch[i][j]->Draw();
      c1->Update();
      //gPad->WaitPrimitive();
      sprintf(name,"h_HGped%d_%d",i,j);
      sprintf(filepath,"../plots/Pedestal/Pedestal_HG_Chip%d_Ch%d.pdf",i,j);
      c1->SetTitle(name);
      //c1->SaveAs(filepath);
              
      h_LGped_ch[i][j]->Draw();
      c1->Update();
      //gPad->WaitPrimitive();
      sprintf(name,"h_LGped%d_%d",i,j);
      sprintf(filepath,"../plots/Pedestal/Pedestal_LG_Chip%d_Ch%d.pdf",i,j);
      c1->SetTitle(name);
      //c1->SaveAs(filepath);
      */
      
       
      for(int k = 0; k < NSCA; ++k){
	avg_HG_SCA  [i][j][k] = (float)sumHG  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_HG_SCA[i][j][k]=sigmaCal(mem_of_SCA[i][j][k],sumHG[i][j][k],sumsqHG[i][j][k]);
	avg_LG_SCA  [i][j][k] = (float)sumLG  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_LG_SCA[i][j][k]=sigmaCal(mem_of_SCA[i][j][k],sumLG[i][j][k],sumsqLG[i][j][k]);

	cout << sigma_HG_SCA[i][j][k] << endl;
	cout << mem_of_SCA[i][j][k] << endl;
	cout << avg_HG_SCA[i][j][k] << endl;

	
	/*
	h_HGped[i][j][k]->Draw();
	c1->Update();
	//gPad->WaitPrimitive();
	sprintf(name,"h_HGped%d_%d_%d",i,j,k);
	sprintf(filepath,"../plots/Pedestal/Pedestal_HG_Chip%d_Ch%d_SCA%d.pdf",i,j,k);
	c1->SetTitle(name);
	c1->SaveAs(filepath);
              
	h_LGped[i][j][k]->Draw();
	c1->Update();
	//gPad->WaitPrimitive();
	sprintf(name,"h_LGped%d_%d_%d",i,j,k);
	sprintf(filepath,"../plots/Pedestal/Pedestal_LG_Chip%d_Ch%d_SCA%d.pdf",i,j,k);
	c1->SetTitle(name);
	c1->SaveAs(filepath);
	*/

	noisy_SCA_check[i][j][k] = (sigma_HG_SCA[i][j][k] - sigma_HG_ch[i][j])/sigma_HG_ch[i][j];
	h_noisy_SCA->Fill(noisy_SCA_check[i][j][k]);
	if(noisy_SCA_check[i][j][k]>1){
	  cout << "chip_" << i << "ch_" << j << "SCA" << k << endl;
	}
	delete h_HGped[i][j][k];
	delete h_LGped[i][j][k];
      }
      if(noisy_ch_check[i][j]>1){
	cout << "chip_" << i << "ch_" << j << endl;
      }
      //outfile->Write();
      delete h_HGped_ch[i][j];
      delete h_LGped_ch[i][j];
    }
  }
  //h_noisy_SCA->Draw();
  //c1->Update();
  //gPad->WaitPrimitive();
  //h_noisy_ch->Draw();
  //c1->Update();
  //gPad->WaitPrimitive();
  
  outfile2->Write();
  outfile2->Close();
  
  delete c1;
  
}

void makePlots::read_P_and_N(string ped_file){
  

  char HG_name[100],LG_name[100];
  sprintf(HG_name,"%s_HG.txt",ped_file.c_str());
  sprintf(LG_name,"%s_LG.txt",ped_file.c_str());
  ifstream inHG(HG_name);
  ifstream inLG(LG_name);
  if( !inHG.is_open() || !inLG.is_open()){
    cout << "File not found! Either" << HG_name << " or " << LG_name
	 << "doesn't exist!" << endl;
    return;}
  else{
    cout << "Input ped file is :" << endl;
    cout << "1. " << HG_name << "\n" << "2. "<< LG_name << endl;
    string line;
    getline(inHG,line); // remove header
    int chip,ch;
    int testeof;
    while(true){
      inHG >> testeof;
      if( inHG.eof() ) break;
      else{
	chip = testeof;
	inHG >> ch;
	for(int sca = 0; sca < NSCA ; ++sca)
	  inHG >> avg_HG_SCA[chip][ch][sca];
	inHG >> chip >> ch;
	for(int sca = 0; sca < NSCA ; ++sca)
	  inHG >> sigma_HG_SCA[chip][ch][sca];
      }      
    }
    getline(inLG,line); // remove header
    while(true){
      inLG >> testeof;
      if( inLG.eof() ) break;
      else{
	chip = testeof;
	inLG >> ch;
	for(int sca = 0; sca < NSCA ; ++sca)
	  inLG >> avg_LG_SCA[chip][ch][sca];
	inLG >> chip >> ch;
	for(int sca = 0; sca < NSCA ; ++sca)
	  inLG >> sigma_LG_SCA[chip][ch][sca];
      }      
    }
    
    cout << "Reading pedestal file done!" << endl;
    inHG.close();
    inLG.close();
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
