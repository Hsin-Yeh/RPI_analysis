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
#include "TStyle.h"

//Constructo
makePlots::makePlots(TChain* inchain):Chain1(inchain)
{
  readmap();
  cout << "Constructor of makePlot ... \n\n" << endl;
}
//Destructor
makePlots::~makePlots()
{
  cout << "\n\n";
  cout << "Destructor of makePlot ... " << endl;
}

void makePlots::Init(){
  yamlReader();
  Chain1->SetBranchAddress("event",&event);
  Chain1->SetBranchAddress("chip",&chip);
  Chain1->SetBranchAddress("roll",&roll);
  Chain1->SetBranchAddress("dacinj",&dacinj);
  Chain1->SetBranchAddress("timesamp",&timesamp);
  Chain1->SetBranchAddress("hg",&hg);
  Chain1->SetBranchAddress("lg",&lg);
  Chain1->SetBranchAddress("tot_fast",&tot_fast);
  Chain1->SetBranchAddress("tot_slow",&tot_slow);
  Chain1->SetBranchAddress("toa_rise",&toa_rise);
  Chain1->SetBranchAddress("toa_fall",&toa_fall);
  //outfile = new TFile("output.root","RECREATE");
  app = new TApplication("app",0,0);
}

void makePlots::Ntuplizer(){


  
}

void makePlots::PlotProducer(){

  
  //==================== Define Parameters ====================
    
  char title[200];
  int TotalEntries = Chain1->GetEntries();
  int Nevents = TotalEntries/NCHIP;
  int cross_num = 6;
  int MaxTS = 3; //choose this time sample to be the peak  
  int dac = -1; int test = 0;
  
  int ADC_H_InjCh[Nevents], ADC_L_InjCh[Nevents], TOTS_InjCh[Nevents]; 
  int ADC_H_Cross[cross_num][Nevents],ADC_L_Cross[cross_num][Nevents];
  int ADC_H_ConnectedCh[NformatCH][Nevents],ADC_L_ConnectedCh[NformatCH][Nevents], TOT_ConnectedCh[NformatCH][Nevents];
  int ADC_H_AllCh[NCHANNEL][Nevents], ADC_L_AllCh[NCHANNEL][Nevents], TOT_AllCh[NCHANNEL][Nevents];
  int ADC_H_NoisyChannel[Nevents];
  int dac_ctrl[Nevents];

  
  // Define Histograms 

  TH1D *h = new TH1D("h","",100,150,250); //("title","",slice,star,end)
  TH1D *h_TOTS = new TH1D("h_TOTS","",100,5,500);
  TH1D *h_TOTF = new TH1D("h_TOTF","",100,1000,3000);
  TH1D *h_TOAR = new TH1D("h_TOAR","",100,1000,3000);
  TH1D *h_TOAF = new TH1D("h_TOAF","",100,1000,3000);
    
  // Set Output Root File
  
  string Input_filename;
  Input_filename = input_RUN;
  int start = Input_filename.find_last_of("/");
  int end   = Input_filename.find(".root");
  string outf = Input_filename.substr(start+1,end-start-1);
  sprintf(title,"output_root/%s.root",outf.c_str());
  TFile *outfile = new TFile(title,"recreate");
  TTree *outT1 = new TTree("Rechit_var","Rechit_var");


  //==================== Initialize ====================
  
  for (int i=0; i<Nevents; i++){
    ADC_H_InjCh[i] = 0;
    ADC_L_InjCh[i] = 0;
    dac_ctrl[i] = i;
    for (int cross_n = 0; cross_n < cross_num; cross_n++){
      ADC_H_Cross[cross_n][i] = 0;
      ADC_L_Cross[cross_n][i] = 0;
    }
  }

  //  yamlReader();
  Crosstalk(Inj_ch);
  
  //==================== Loop over the events ====================
   
  for(int entry = 0; entry < TotalEntries ; ++entry){  
    Chain1 -> GetEntry(entry);
    dac_ctrl[event] = dacinj;
    
    // Filling hg, lg data (with 13 SCA)
    
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = timesamp[i];    
    for(int sca=0; sca<NSCA; sca++){
      if(TS[sca]==MaxTS){
	
	ADC_H_InjCh[event] = hg[sca][Inj_ch]; // Filling the Inj_ch
	ADC_L_InjCh[event] = lg[sca][Inj_ch];
	
	for(int ch = 0; ch < 32; ch++){
	    ADC_H_ConnectedCh[ch+chip*32][event] = hg[sca][ch*2]; // Filling all the connected channels
	    ADC_L_ConnectedCh[ch+chip*32][event] = lg[sca][ch*2];
	    cout << hg[sca][ch*2] << endl;
	}
	
	for(int icross = 0; icross < cross_num; icross++){
	  ADC_H_Cross[icross][event] = hg[sca][cross_ch[icross]]; // Filling FirstRing around Inj_ch
	  ADC_L_Cross[icross][event] = lg[sca][cross_ch[icross]];
	}
	
      }
    }

    // Filling tot, toa data (without SCA)
    for(int ch = 0; ch < NformatCH/NCHIP; ch++){
      TOT_ConnectedCh[ch+chip*32][event] = tot_slow[ch*2];
    }
    
  }

  
  //==================== End of Loop ====================

  
  
  //==================== Fit & Plots ====================
  

  // Define Plotting Paramter
  
  char plot_title[200], leg[50], img_title[50];
  sprintf(title,"plots/TBHexaboard/module%d",ModuleNumber);
  string plotfolder_path(title);

  TLegend *legend = new TLegend(0.85,0.8,1.,1.);
  //legend->SetNColumns(2);
  TCanvas* c1 = new TCanvas();
  TCanvas* c2 = new TCanvas("c2","c2",6400,3600);
  c2->Divide(8,8);
  gStyle->SetOptStat(0);


  // Define Fitting Parameter

  double slope_h[NformatCH], slope_l[NformatCH], slope_tot[NformatCH];
  double slope_h_InjCh, slope_l_InjCh;


  // Define TGraphs

  TGraph* gh = new TGraph(Nevents,dac_ctrl,ADC_H_InjCh);
  TGraph* gl = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh);
  TGraph* gTOT = new TGraph(Nevents,dac_ctrl,TOTS_InjCh);
  TGraph* ghlratio = new TGraph(Nevents,ADC_L_InjCh,ADC_H_InjCh);
  TGraph* gltratio = new TGraph(Nevents,TOTS_InjCh,ADC_L_InjCh);
  TGraph** gcross_h = new TGraph*[cross_num];
  TGraph** gcross_l = new TGraph*[cross_num];
  TGraph** gcross_tot = new TGraph*[cross_num];
  TGraph** gh_ConnectedCh = new TGraph*[NformatCH];
  TGraph** gl_ConnectedCh = new TGraph*[NformatCH];
  TGraph** gtot_ConnectedCh = new TGraph*[NformatCH];
  TGraph* gnoisy_h = new TGraph(Nevents,dac_ctrl,ADC_H_NoisyChannel);
  TGraph* gcorrelation_l = new TGraph(Nevents,ADC_L_InjCh,ADC_H_Cross[0]);
  TMultiGraph* multig_cross_h = new TMultiGraph();
  TMultiGraph* multig_cross_l = new TMultiGraph();


  // Plot !!!!

  for(int ch = 0; ch < NformatCH; ch++){
    c2->cd(ch+1);
    gh_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_ConnectedCh[ch]);
    sprintf(plot_title,"HG %d",ch);
    gh_ConnectedCh[ch]->SetTitle(plot_title);
    gh_ConnectedCh[ch]->GetXaxis()->SetTitle("DAC");
    gh_ConnectedCh[ch]->GetYaxis()->SetTitle("ADC");
    gh_ConnectedCh[ch]->SetMarkerStyle(7);
    gh_ConnectedCh[ch]->Fit("pol1","","",1000,4000);
    //gh_ConnectedCh[ch]->Fit("pol1");
    //    gh_ConnectedCh[ch]->Draw("AP");
    
    gl_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_ConnectedCh[ch]);
    sprintf(plot_title,"LG ch = %d",ch);
    gl_ConnectedCh[ch]->SetTitle(plot_title);
    //gl_ConnectedCh[ch]->GetYaxis()->SetRangeUser(0,1500);
    gl_ConnectedCh[ch]->GetXaxis()->SetTitle("DAC");
    gl_ConnectedCh[ch]->GetYaxis()->SetTitle("ADC");
    gl_ConnectedCh[ch]->SetMarkerStyle(7);
    gl_ConnectedCh[ch]->Fit("pol1","","",1000,4000);
    //gl_ConnectedCh[ch]->Fit("pol1");
    gl_ConnectedCh[ch]->Draw("AL");

    gtot_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_ConnectedCh[ch]);
    sprintf(plot_title,"tot ch = %d",ch);
    gtot_ConnectedCh[ch]->SetTitle(plot_title);
    //gtot_ConnectedCh[ch]->GetYaxis()->SetRangeUser(0,1500);
    gtot_ConnectedCh[ch]->GetXaxis()->SetTitle("DAC");
    gtot_ConnectedCh[ch]->GetYaxis()->SetTitle("ADC");
    gtot_ConnectedCh[ch]->SetMarkerStyle(7);
    gtot_ConnectedCh[ch]->Fit("pol1","","",1000,4000);
    //gtot_ConnectedCh[ch]->Fit("pol1");
    gtot_ConnectedCh[ch]->Draw("AL");
    
    
    TF1* Linear_fit_h = gh_ConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l = gl_ConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot = gtot_ConnectedCh[ch]->GetFunction("pol1");
    slope_h[ch] = Linear_fit_h->GetParameter(1);
    slope_l[ch] = Linear_fit_l->GetParameter(1);
    slope_tot[ch] = Linear_fit_tot->GetParameter(1);
  }

  //c2->Update();
  //gPad->WaitPrimitive();

  delete c2;

  c1->cd();
  //gh_ConnectedCh[18]->Draw("AP");
  //c1->Update();
  //gPad->WaitPrimitive();

  gh->Fit("pol1","","",0,200);
  TF1* Linear_fit_h_InjCh = gh->GetFunction("pol1");
  slope_h_InjCh = Linear_fit_h_InjCh->GetParameter(1);
  sprintf(plot_title,"Inj_CH_Chip0Channel%dTS%d_HG",Inj_ch,MaxTS);
  gh->SetTitle(plot_title);
  gh->GetXaxis()->SetTitle("DAC");
  gh->GetYaxis()->SetTitle("ADC");
  gh->SetMarkerStyle(7);
  gh->Draw("AP");
  c1->Update();
  sprintf(title,"%s/%s.pdf",plotfolder_path.c_str(),plot_title);
  //c1->SaveAs(title);
  //gPad->WaitPrimitive();

  
  gl->Fit("pol1","","",0,500);
  TF1* Linear_fit_l_InjCh = gl->GetFunction("pol1");
  slope_l_InjCh = Linear_fit_l_InjCh->GetParameter(1);
  sprintf(plot_title,"Inj_CH_Chip0Channel%dTS%d_LG",Inj_ch,MaxTS);
  gl->SetTitle(plot_title);
  gl->GetXaxis()->SetTitle("DAC");
  gl->GetYaxis()->SetTitle("ADC");
  gl->SetMarkerStyle(7);
  gl->Draw("AP");
  c1->Update();
  sprintf(title,"%s/%s.pdf",plotfolder_path.c_str(),plot_title);
  //c1->SaveAs(title);
  //gPad->WaitPrimitive();

  sprintf(plot_title,"Low Gain correlation");
  gcorrelation_l->SetTitle(plot_title);
  gcorrelation_l->GetXaxis()->SetTitle("CH2");
  gcorrelation_l->GetYaxis()->SetTitle("CH0");
  gcorrelation_l->SetMarkerStyle(7);
  gcorrelation_l->Draw("AP");
  c1->Update();
  //gPad->WaitPrimitive();

  for(int cross_n = 0; cross_n < cross_num; cross_n++){
    gcross_h[cross_n] = new TGraph(Nevents,dac_ctrl,ADC_H_Cross[cross_n]);
    gcross_l[cross_n] = new TGraph(Nevents,dac_ctrl,ADC_L_Cross[cross_n]);
      
    sprintf(plot_title,"CH %d High Gain",cross_ch[cross_n]);
    gcross_h[cross_n]->SetTitle(plot_title);
    gcross_h[cross_n]->GetXaxis()->SetTitle("DAC");
    gcross_h[cross_n]->GetYaxis()->SetTitle("ADC");
    gcross_h[cross_n]->SetMarkerStyle(26);
    gcross_h[cross_n]->SetMarkerSize(0.4);
    gcross_h[cross_n]->SetMarkerColor(cross_n);
    gcross_h[cross_n]->Draw("AP");
    c1->Update();
    //gPad->WaitPrimitive();
      
    sprintf(plot_title,"CH %d Low Gain",cross_ch[cross_n]);
    gcross_l[cross_n]->SetTitle(plot_title);
    gcross_l[cross_n]->GetXaxis()->SetTitle("DAC");
    gcross_l[cross_n]->GetYaxis()->SetTitle("ADC");
    gcross_l[cross_n]->SetMarkerSize(0.4);
    gcross_l[cross_n]->SetMarkerStyle(26);
    gcross_l[cross_n]->Draw("AP");      
    c1->Update();
    //gPad->WaitPrimitive();

    multig_cross_h->Add(gcross_h[cross_n]);
    multig_cross_l->Add(gcross_l[cross_n]);
  }  
  sprintf(plot_title,"FirstRing_&_UnConnectedCh_around_Chip0Channel%dTS%dHG",Inj_ch,MaxTS);
  multig_cross_h->SetTitle(plot_title);
  multig_cross_h->Draw("AP");
  multig_cross_h->GetXaxis()->SetTitle("DAC");
  multig_cross_h->GetYaxis()->SetTitle("ADC");
  c1->BuildLegend(0.1,0.7,0.3,0.9);
  c1->Update();
  sprintf(title,"%s/%s.pdf",plotfolder_path.c_str(),plot_title);
  //c1->SaveAs(title);
  //gPad->WaitPrimitive();
  

  //************************************************** Cross talk Plots **************************************************//

  int NNoisyCh = 7;
  int NoisyChannel[7] = {248,186,214,120,126,42,254};
  
  TH2Poly *polyh = new TH2Poly;
  TH2Poly *polyl = new TH2Poly;
  InitTH2Poly(*polyh);
  InitTH2Poly(*polyl);
  
  for(int ch = 0; ch < NformatCH; ch++){
    float X, Y;
    bool NoisyBool = false;
    X = CHmap[ch].first;
    Y = CHmap[ch].second;
    if(ch%32==Inj_ch/2){
      polyh->Fill(X,Y,0.5);
      polyl->Fill(X,Y,0.05);
    }
    else {
      for(int iNoisy = 0; iNoisy < NNoisyCh; iNoisy++){
	if(ch == NoisyChannel[iNoisy]/2) {NoisyBool = true;}
      }
      if(!NoisyBool){
	polyh->Fill(X,Y,slope_h[ch]);
	polyl->Fill(X,Y,slope_l[ch]);
      }
    }
  }
  
  sprintf(plot_title,"Slope of HG TS%d vs Injdac, Inj_ch=%d",MaxTS,Inj_ch);
  polyh->SetTitle(plot_title);
  polyh->Draw("colztext");
  c1->Update();
  sprintf(title,"%s/%s.pdf",plotfolder_path.c_str(),plot_title);
  c1->SaveAs(title);
  //  gPad->WaitPrimitive();

  sprintf(plot_title,"Slope of LG TS%d vs Injdac, Inj_ch=%d",MaxTS,Inj_ch);
  polyl->SetTitle(plot_title);
  polyl->Draw("colztext");
  c1->Update();
  sprintf(title,"%s/%s.pdf",plotfolder_path.c_str(),plot_title);
  c1->SaveAs(title);
  //gPad->WaitPrimitive();

  /*
  InitTH2Poly(*polyh);
p  InitTH2Poly(*polyl);
  
  for(int ch = 0; ch < NformatCH; ch++){
    float X, Y;
    X = CHmap[ch].first;
    Y = CHmap[ch].second;    
    if(ch%32==Inj_ch/2){
      polyh->Fill(X,Y,0.8);
      polyl->Fill(X,Y,0.2);
    }
    else {
      polyh->Fill(X,Y,slope_h[ch]);
      polyl->Fill(X,Y,slope_l[ch]);
    }
  }
  */
  
  delete c1;  
}


//****************************************************************************************************//
//****************************************************************************************************//

/*
void makePlots::Evt_display(){

  //  Init();
  TCanvas* c1 = new TCanvas();
  int MaxTS = 5;

  int Nevents = Chain1->GetEntries();
  char plot_title[50];
  for(int ev = 0; ev < Nevents ; ++ev){
    if(ev % 10 != 0) continue;

    Chain1 -> GetEntry(ev); // Get HITcollection from root file event
    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);
    
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITCOLLECTION->rollposition[i];

    int nhits = HITCOLLECTION->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){  
      H = HITCOLLECTION->Hits.at(hit);
      if(!H.CCorNC) continue; // remove unconnected channel
      int forCH = H.chip*32+H.ch/2;
      float X = CHmap[forCH].first;
      float Y = CHmap[forCH].second;	  

      for(int sca = 0; sca < NSCA; ++sca){
p	if(TS[sca] == MaxTS ){
	  //cout << H.SCA_lg[sca] << endl;
	  //	  H.SCA_hg[sca] -= avg_HG[H.chip][H.ch][sca]; // pedestal subtraction
	  //poly->Fill(X,Y,H.SCA_lg[sca]);
	  poly->Fill(X,Y,H.SCA_lg[sca]);
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
  
  delete c1;
}
*/


//****************************************************************************************************//
//****************************************************************************************************//



void makePlots::yamlReader(){

  int start, end;
  string rootFileName(input_RUN);
  start = rootFileName.find("/unpack");
  string temp = rootFileName.replace(start+1,13,"yaml");
  start = temp.find("_pedestal.root");
  string yamlFileName = temp.replace(start,14,".yaml");
  
  string searchstr;
  ifstream yamlFile(yamlFileName);
  if(!yamlFile.is_open()){
    cout << "Did not find injection file " << yamlFileName
	 << ".\n Take this run as pedestal.(Inj_dac = 0)" << endl;
  }
  if(yamlFile.is_open()){
    for(int header = 0; header < 3; header++) {getline(yamlFile,searchstr);}
    start = searchstr.find("[");
    searchstr = searchstr.substr(start+1,start+2);
    end = searchstr.find("]");
    searchstr = searchstr.erase(end);
    Inj_ch = atoi(searchstr.c_str());
    cout << Inj_ch << endl;

    for(int header = 0; header < 13; header++) {getline(yamlFile,searchstr);}
    start = searchstr.find("'");
    searchstr = searchstr.substr(start+1,start+2);
    end = searchstr.find("',");
    searchstr = searchstr.erase(end);
    ModuleNumber = atoi(searchstr.c_str());
    cout << ModuleNumber << endl;

    
  }

}


//****************************************************************************************************//
//****************************************************************************************************//



void makePlots::Inj_Pulse_display(){

  //Init();
  TGraph *gr;
  int Nevents = Chain1->GetEntries();
  char plot_title[50];
  TCanvas* c1 = new TCanvas();
  int lg_transpose[64][13];
  
  for(int ev = 0; ev < Nevents ; ++ev){
    if(ev % 1 != 0) continue;
    Chain1 -> GetEntry(ev);
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = timesamp[i];
    
    for (int sca = 0; sca < NSCA; sca++){
      for(int ch = 0; ch < 64; ch++){
	lg_transpose[ch][sca] = lg[sca][ch];
      }
    }
    
    
    gr = new TGraph(NSCA, TS, lg_transpose[Inj_ch] );
    gr->SetMarkerColor(chip+2);
    gr->SetMarkerStyle(22);
    gr->SetMarkerSize(1.2);
    gr->Draw("AP");
    sprintf(plot_title,"LG_evt %d chip %d",ev,chip);
    gr->SetTitle(plot_title);
    gr->GetXaxis()->SetTitle("TS");
    gr->GetYaxis()->SetTitle("ADC");

    c1->Update();
    //if(ev == 400){
    //sprintf(plot_title,"%s.png",plot_title);
    //c1->SaveAs(plot_title);} // remove the comment to save plots
    gPad->WaitPrimitive();
      
  }
  delete c1;
  
}


//****************************************************************************************************//
//****************************************************************************************************//


  
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


//****************************************************************************************************//
//****************************************************************************************************//


void makePlots::Crosstalk(Int_t CH){

  TCanvas* c1 = new TCanvas();
  int cross_num = 6;
  int formatInj_Ch = CH/2;
  float Xdist = 0.974452;
  float Ydist = 0.5626;
  float cross_posx[cross_num], cross_posy[cross_num];
  float X = CHmap[formatInj_Ch].first;
  float Y = CHmap[formatInj_Ch].second;
  
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
      //  cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
      ch++;
      if(ch>256) break;
    }
    cross_ch[i] = (ch-chip*32)*2;
  }
  /*
    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly);
  
    for(int i=0; i<6; i++){
    cout << cross_ch[i] << endl;
    forCH = chip*32+cross_ch[i]/2;
    X = CHmap[forCH].first;
    Y = CHmap[forCH].second;	      
    poly->Fill(X,Y,cross_ch[i]);
    poly->Draw("colztext0");
    c1->Update();
    //gPad->WaitPrimitive();
    }
  */
  delete c1;
  
}


//****************************************************************************************************//
//****************************************************************************************************//

/*
void makePlots::Crosstalk_2ndRing(Int_t CH){

  TCanvas* c1 = new TCanvas();
  int chip = 0;
  int cross_num_2ndRing = 12;
  int forCH = chip*32+CH/2;
  float Xdist = 0.974452;
  float Ydist = 0.5626;
  float cross_posx[cross_num_2ndRing], cross_posy[cross_num_2ndRing];
  float X = CHmap[forCH].first;
  float Y = CHmap[forCH].second;
  
  cross_posx[0] = X + 2*Xdist;
  cross_posy[0] = Y;
  cross_posx[1] = X + 2*Xdist;
  cross_posy[1] = Y - 2*Ydist;
  cross_posx[2] = X + Xdist;
  cross_posy[2] = Y - 3*Ydist;
  cross_posx[3] = X;
  cross_posy[3] = Y - 4*Ydist;
  cross_posx[4] = X - Xdist;
  cross_posy[4] = Y - 3*Ydist;
  cross_posx[5] = X - 2*Xdist;
  cross_posy[5] = Y - 2*Ydist;
  cross_posx[6] = X - 2*Xdist;
  cross_posy[6] = Y;
  cross_posx[7] = X - 2*Xdist;
  cross_posy[7] = Y + 2*Ydist;
  cross_posx[8] = X - Xdist;
  cross_posy[8] = Y + 3*Ydist;
  cross_posx[9] = X;
  cross_posy[9] = Y + 4*Ydist;
  cross_posx[10] = X + Xdist;
  cross_posy[10] = Y + 3*Ydist;
  cross_posx[11] = X + 2*Xdist;
  cross_posy[11] = Y + 2*Ydist;


  for(int i=0; i<cross_num_2ndRing; i++){
    //cout << cross_posx[i] << " " << cross_posy[i] << endl;
    int ch = 0;
    
    while(abs(CHmap[ch].first-cross_posx[i]) > 1e-4 || abs(CHmap[ch].second-cross_posy[i]) > 1e-4){
      //  cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
      ch++;
      if(ch>256) break;
    }
    cross_ch_2ndRing[i] = (ch-chip*32)*2;
    //    cross_ch[i] = ch;
  }

  TH2Poly *poly = new TH2Poly;
  InitTH2Poly(*poly);

  for(int i=0; i<6; i++){
    cout << cross_ch_2ndRing[i] << endl;
    forCH = chip*32+cross_ch_2ndRing[i]/2;
    X = CHmap[forCH].first;
    Y = CHmap[forCH].second;	      
    poly->Fill(X,Y,forCH);
    poly->Draw("colztext0");
    c1->Update();
    //gPad->WaitPrimitive();
  }
  
  delete c1;
  
}
*/

//****************************************************************************************************//
//****************************************************************************************************//


void makePlots::readmap(){
  ifstream file("./src_txtfile/CH_map.txt");
  string line;
  int ichip,ich,itype,iformatCH;
  double iposx, iposy;
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> ichip >> ich >> iposx >> iposy >> itype;
    iformatCH = ichip*32 + ich/2;
    CHmap[iformatCH] = make_pair(iposx,iposy);}
  file.close();
  //Since there is no such pad, assign a unreasonable value
  CHmap[2*32+60/2] = make_pair(1000.,1000.);

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
    poly.AddBin(CellXYsize, HexX, HexY);
  }
  file.close();
}
