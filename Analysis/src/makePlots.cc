#include "makePlots.h"
#include "PlotSetting.h"
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
#include "TSystem.h"
#include "TImage.h"
#include "TStyle.h"
#include "TExec.h"

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



void makePlots::Gain_factor_producer(){
  
  //-------------------- Define Parameters --------------------

  int TotalEntries = Chain1->GetEntries();
  int Nevents = TotalEntries/NCHIP;
  int MaxTS = 2;

  int ADC_H_InjCh_Chip[NCHIP][Nevents], ADC_L_InjCh_Chip[NCHIP][Nevents], TOT_InjCh_Chip[NCHIP][Nevents];
  int dac_ctrl[Nevents];
  int HGLGfitmax[NCHIP], TOTLGfitmax[NCHIP];
  bool HGTP_flag[NCHIP], LGTP_flag[NCHIP];

  //==================== Set Output File ====================

  


  //==================== Initialize ====================

  for(int ichip = 0; ichip < NCHIP; ichip++){
    HGTP_flag[ichip] = false;
    LGTP_flag[ichip] = false;
  }


  //==================== Loop Over Events ====================

  for(int ev = 0; ev < TotalEntries; ev++){
    if(ev%1000 == 0){ cout << "Now Processing = " << ev << endl;}
    Chain1->GetEntry(ev);
    dac_ctrl[event] = dacinj;
    
    int TS0_sca, MaxTS_sca;
    for(int sca = 0 ; sca < NSCA ; sca++) {
      TS[sca] = timesamp[sca];
      if (timesamp[sca] == 0) { TS0_sca = sca ; }
      if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
    }
    
    ADC_H_InjCh_Chip[chip][event] = ( hg[MaxTS_sca][Inj_ch] - hg[TS0_sca][Inj_ch] );
    ADC_L_InjCh_Chip[chip][event] = ( lg[MaxTS_sca][Inj_ch] - lg[TS0_sca][Inj_ch] );
    TOT_InjCh_Chip[chip][event] = tot_slow[Inj_ch];

    if(ADC_H_InjCh_Chip[chip][event] > HGTP && HGTP_flag[chip] == false){
      HGTP_flag[chip] = true;
      HGLGfitmax[chip] = ADC_L_InjCh_Chip[chip][event];
    }
    if(ADC_L_InjCh_Chip[chip][event] > LGTP && LGTP_flag[chip] == false){
      LGTP_flag[chip] = true;
      TOTLGfitmax[chip] = TOT_InjCh_Chip[chip][event];
    }
  }

  
    
  //==================== End Loop ====================

  //...

  //==================== Plots ====================

  char pltTit[100];
  string Xtit, Ytit, Opt;
  int MkSty, MkClr, MkSize, LClr, LWid, fitmin, fitmax;
  bool Stat, Wait, SavePlot;

  TGraph** gl = new TGraph*[NCHIP];
  TGraph** gtot = new TGraph*[NCHIP];
  TGraph** LG2HG = new TGraph*[NCHIP];
  TGraph** TOT2LG = new TGraph*[NCHIP];

  for(int ichip = 0; ichip < NCHIP; ichip++){
    gl[ichip] = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh_Chip[ichip]);
    gtot[ichip] = new TGraph(Nevents,dac_ctrl,TOT_InjCh_Chip[ichip]);
    LG2HG[ichip] = new TGraph(Nevents,ADC_L_InjCh_Chip[ichip],ADC_H_InjCh_Chip[ichip]);
    TOT2LG[ichip] = new TGraph(Nevents,TOT_InjCh_Chip[ichip],ADC_L_InjCh_Chip[ichip]);
    
    LG2HG[ichip]->Fit("pol1","","",fitmin = 0,fitmax = HGLGfitmax[ichip]);
    TOT2LG[ichip]->Fit("pol1","","",fitmin = TOTOffSet,fitmax = TOTLGfitmax[ichip]);
    
    TF1* Linear_fit_LG2HG = LG2HG[ichip]->GetFunction("pol1");
    TF1* Linear_fit_TOT2LG = TOT2LG[ichip]->GetFunction("pol1");
    LG2HG_Conversion[ichip] = Linear_fit_LG2HG->GetParameter(1);
    TOT2LG_Conversion[ichip] = Linear_fit_TOT2LG->GetParameter(1);

    sprintf(pltTit,"LG_Chip%d",ichip);
    Plot.TGraphPlotSetting(*gl[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 25, MkClr = 1+ichip, MkSize = 0.5, LClr = 1+ichip, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    sprintf(pltTit,"TOT_Chip%d",ichip);
    Plot.TGraphPlotSetting(*gtot[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 25, MkClr = 1+ichip, MkSize = 0.5, LClr = 1+ichip, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    
    sprintf(pltTit,"LG2HG_Chip%d",ichip);
    Plot.TGraphPlotSetting(*LG2HG[ichip], pltTit, Xtit = "LG", Ytit = "HG",
		      MkSty = 25, MkClr = 1+ichip, MkSize = 0.5, LClr = 1+ichip, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    
    sprintf(pltTit,"TOT2LG_Chip%d",ichip);
    Plot.TGraphPlotSetting(*TOT2LG[ichip], pltTit, Xtit = "TOT", Ytit = "LG",
		      MkSty = 25, MkClr = 1+ichip, MkSize = 0.5, LClr = 1+ichip, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
  }
}



//****************************************************************************************************//
//                                        PlotProducer                                                //
//****************************************************************************************************//




void makePlots::PlotProducer(){

  
  //==================== Define Parameters ====================
     
  char title[200];
  int TotalEntries = Chain1->GetEntries();
  int Nevents = TotalEntries/NCHIP;
  int cross_num = 6;
  int MaxTS = 2; //choose this time sample to be the peak  
  int dac = -1; int test = 0;
  
  int ADC_H_InjCh[Nevents], ADC_L_InjCh[Nevents], TOTS_InjCh[Nevents];
  int ADC_H_InjCh_Chip[NCHIP][Nevents], ADC_L_InjCh_Chip[NCHIP][Nevents], TOT_InjCh_Chip[NCHIP][Nevents];
  int ADC_H_Cross[cross_num][Nevents],ADC_L_Cross[cross_num][Nevents], ADC_H_Cross_Chip[NCHIP][cross_num][Nevents], ADC_L_Cross_Chip[NCHIP][cross_num][Nevents];
  int ADC_H_FirstRing[NCHIP][Nevents],ADC_L_FirstRing[NCHIP][Nevents];
  int ADC_H_ConnectedCh[NformatCH][Nevents],ADC_L_ConnectedCh[NformatCH][Nevents], TOT_ConnectedCh[NformatCH][Nevents];
  int ADC_H_AllCh[NCHANNEL][Nevents], ADC_L_AllCh[NCHANNEL][Nevents], TOT_AllCh[NCHANNEL][Nevents];
  int ADC_H_NoisyChannel[Nevents];
  double ADC_L_InjCh_Chip_double[NCHIP][Nevents], Ratio_L_FirstRing_InjCh[NCHIP][event], ADC_H_InjCh_Chip_double[NCHIP][Nevents], Ratio_H_FirstRing_InjCh[NCHIP][event];
  int dac_ctrl[Nevents];
  int cross_ch_chip[cross_num][NCHIP];
  bool cross_type_chip[cross_num][NCHIP];


  
  
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
    for(int ichip = 0; ichip < NCHIP; ichip++){
      ADC_L_FirstRing[ichip][i] = 0;
    }
    for (int cross_n = 0; cross_n < cross_num; cross_n++){
      ADC_H_Cross[cross_n][i] = 0;
      ADC_L_Cross[cross_n][i] = 0;
    }
  }

  for(int ichip = 0; ichip < NCHIP; ichip++){
    Crosstalk(ichip,Inj_ch);
    for(int cross_n = 0; cross_n < cross_num; cross_n++){
      cross_ch_chip[cross_n][ichip] = cross_ch[cross_n];
      cross_type_chip[cross_n][ichip] = cross_type[cross_n];
    }
  }
    
  
  //==================== Loop over the events ====================
   
  for(int entry = 0; entry < TotalEntries ; ++entry){
    
    if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }    
    Chain1 -> GetEntry(entry);
    dac_ctrl[event] = dacinj;
    
    // Filling hg, lg data (with 13 SCA)
    int TS0_sca, MaxTS_sca;
    for(int sca = 0 ; sca < NSCA ; sca++) {
      TS[sca] = timesamp[sca];
      if (timesamp[sca] == 0) { TS0_sca = sca ; }
      if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
    }
	
    ADC_H_InjCh[event] = hg[MaxTS_sca][Inj_ch]; // Filling the Inj_ch
    ADC_L_InjCh[event] = lg[MaxTS_sca][Inj_ch];

    ADC_H_InjCh_Chip[chip][event] = hg[MaxTS_sca][Inj_ch];
    ADC_L_InjCh_Chip[chip][event] = lg[MaxTS_sca][Inj_ch];
    ADC_H_InjCh_Chip_double[chip][event] = ( hg[MaxTS_sca][Inj_ch] - hg[TS0_sca][Inj_ch] );
    ADC_L_InjCh_Chip_double[chip][event] = ( lg[MaxTS_sca][Inj_ch] - lg[TS0_sca][Inj_ch] );
	
    for(int ch = 0; ch < 32; ch++){
      ADC_H_ConnectedCh[ch+chip*32][event] = hg[MaxTS_sca][ch*2]; // Filling all the connected channels
      ADC_L_ConnectedCh[ch+chip*32][event] = lg[MaxTS_sca][ch*2];
      ADC_H_AllCh[ch+chip*64][event] = hg[MaxTS_sca][ch];
      ADC_H_AllCh[ch+chip*64][event] = hg[MaxTS_sca][ch];
    }
	
    for(int icross = 0; icross < cross_num; icross++){
      if(cross_type_chip[icross][chip] == false) break;
      ADC_H_Cross[icross][event] = hg[MaxTS_sca][cross_ch_chip[chip][icross]]; // Filling FirstRing around Inj_ch
      ADC_L_Cross[icross][event] = lg[MaxTS_sca][cross_ch_chip[chip][icross]];
      ADC_H_Cross_Chip[chip][icross][event] = hg[MaxTS_sca][cross_ch_chip[chip][icross]];
      ADC_L_Cross_Chip[chip][icross][event] = lg[MaxTS_sca][cross_ch_chip[chip][icross]];
      ADC_H_FirstRing[chip][event] += (hg[MaxTS_sca][cross_ch_chip[chip][icross]] - hg[TS0_sca][cross_ch_chip[chip][icross]]);
      ADC_L_FirstRing[chip][event] += (lg[MaxTS_sca][cross_ch_chip[chip][icross]] - lg[TS0_sca][cross_ch_chip[chip][icross]]);
    }
    Ratio_L_FirstRing_InjCh[chip][event] = ADC_L_FirstRing[chip][event]/(ADC_L_InjCh_Chip_double[chip][event]+ADC_L_FirstRing[chip][event]);
    Ratio_H_FirstRing_InjCh[chip][event] = ADC_H_FirstRing[chip][event]/(ADC_H_InjCh_Chip_double[chip][event]+ADC_H_FirstRing[chip][event]);
	

    // Filling tot, toa data (without SCA)
    for(int ch = 0; ch < NformatCH/NCHIP; ch++){
      TOT_ConnectedCh[ch+chip*32][event] = tot_slow[ch*2];
    }

    TOT_InjCh_Chip[chip][event] = tot_slow[Inj_ch];
    
  }
  
  //... ==================== End of Loop ==================== ...

  
  
  //... ==================== Fit & Plots ==================== ...
  

  // Define Plotting Paramter
  
  char pltTit[200], leg[50], img_title[50];
  string Xtit, Ytit, Opt;
  int MkSty, MkClr, MkSize, LClr, LWid, fitmin, fitmax;
  bool Stat, Wait, SavePlot;


  // Define Fitting Parameter

  double slope_h[NformatCH], slope_l[NformatCH], slope_tot[NformatCH];
  double slope_h_Uncnct[NformatCH], slope_l_Uncnct[NformatCH], slope_tot_Uncnct[NformatCH];
  double slope_h_InjCh, slope_l_InjCh, slope_h_chip[NCHIP], slope_l_chip[NCHIP];
  double CnctID[NformatCH], UncnctID[NformatCH];
  

  for(int ch = 0; ch < NformatCH; ch++){
    UncnctID[ch] = ch*2 + 1;
  }


  // Define TGraphs

  TGraph* gh = new TGraph(Nevents,dac_ctrl,ADC_H_InjCh);
  TGraph* gl = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh);
  TGraph* gTOT = new TGraph(Nevents,dac_ctrl,TOTS_InjCh);
  TGraph* ghlratio = new TGraph(Nevents,ADC_L_InjCh,ADC_H_InjCh);
  TGraph* gltratio = new TGraph(Nevents,TOTS_InjCh,ADC_L_InjCh);
  TGraph* gh_Uncnct_Slope;
  TGraph* gh_Cnct_Slope;
  TGraph** gh_chip = new TGraph*[NCHIP];
  TGraph** gl_chip = new TGraph*[NCHIP];
  TGraph** gtot_chip = new TGraph*[NCHIP];
  TGraph** gh_UnConnectedCh = new TGraph*[NformatCH];
  TGraph** gl_UnConnectedCh = new TGraph*[NformatCH];
  TGraph** gtot_UnConnectedCh = new TGraph*[NformatCH];  
  TGraph** gcross_h = new TGraph*[cross_num];
  TGraph** gcross_l = new TGraph*[cross_num];
  TGraph** gcross_tot = new TGraph*[cross_num];
  TGraph** gh_ConnectedCh = new TGraph*[NformatCH];
  TGraph** gl_ConnectedCh = new TGraph*[NformatCH];
  TGraph** gtot_ConnectedCh = new TGraph*[NformatCH];
  TGraph* gnoisy_h = new TGraph(Nevents,dac_ctrl,ADC_H_NoisyChannel);
  TGraph* gcorrelation_l = new TGraph(Nevents,ADC_L_InjCh,ADC_H_Cross[0]);
  TGraph** gratioRing1_Injch_l = new TGraph*[NCHIP];
  TMultiGraph* multig_cross_h = new TMultiGraph();
  TMultiGraph* multig_cross_l = new TMultiGraph();
  TMultiGraph* multig_InjCh_Chip_h = new TMultiGraph();
  TMultiGraph* multig_InjCh_Chip_l = new TMultiGraph();
  TMultiGraph* multig_InjCh_Chip_hltot = new TMultiGraph();
  TMultiGraph* multig_AllCh_Slope = new TMultiGraph();


  // Plots!!!!!


  // ------------------------------ Fit ------------------------------ //

  // UnconnectedCh
  
  for(int ch = 0; ch < NformatCH; ch++){    

    fitmin = 1300;
    fitmax = 3300;    

    gh_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2+1]);
    gh_UnConnectedCh[ch]->Fit("pol1","","",1300,3300);
    sprintf(pltTit,"HG_%d",ch);
    Plot.TGraphPlotSetting(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2+1]);
    gl_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.TGraphPlotSetting(*gl_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2+1]);
    gtot_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.TGraphPlotSetting(*gtot_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    TF1* Linear_fit_h_Uncnct = gh_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l_Uncnct = gl_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot_Uncnct = gtot_UnConnectedCh[ch]->GetFunction("pol1");
    slope_h_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_l_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_tot_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);

    if(ch%50==0){
      sprintf(pltTit,"NcnctCh%d_InjCh%d_fit%d-%d",ch*2+1,Inj_ch,fitmin,fitmax);
      Plot.TGraphPlotSetting(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
			MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    }
  }
  
  //...Slope vs Uncnct Channel
  
  gh_Uncnct_Slope = new TGraph(NformatCH,UncnctID,slope_h_Uncnct);
  sprintf(pltTit,"UnconnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
  Plot.TGraphPlotSetting(*gh_Uncnct_Slope, pltTit, Xtit = "CH", Ytit = "P1",
		    MkSty = 22, MkClr = 2, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);


  //... ConnectedCh

  for(int ch = 0; ch < NformatCH; ch++){
    fitmin = 1000;
    fitmax = 4000;
    
    gh_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2]);
    gh_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"HG_%d",ch);
    Plot.TGraphPlotSetting(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2]);
    gl_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.TGraphPlotSetting(*gl_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2]);
    gtot_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.TGraphPlotSetting(*gtot_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    TF1* Linear_fit_h = gh_ConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l = gl_ConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot = gtot_ConnectedCh[ch]->GetFunction("pol1");
    slope_h[ch] = Linear_fit_h->GetParameter(1);
    slope_l[ch] = Linear_fit_l->GetParameter(1);
    slope_tot[ch] = Linear_fit_tot->GetParameter(1);

    CnctID[ch] = ch*2;

    if(ch%50==0){
      sprintf(pltTit,"NcnctCh%d_InjCh%d_fit%d-%d",ch*2+1,Inj_ch,fitmin,fitmax);
      Plot.TGraphPlotSetting(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
			MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    }
  }

  //...Slope vs cnct Channel
  
  gh_Cnct_Slope = new TGraph(NformatCH,CnctID,slope_h);
  sprintf(pltTit,"ConnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
  Plot.TGraphPlotSetting(*gh_Cnct_Slope, pltTit, Xtit = "Ch", Ytit = "P1",
		    MkSty = 22, MkClr = 3, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);


  //------------------------------ Change fitting Range ------------------------------//

  //... UnconncetedCh

  for(int ch = 0; ch < NformatCH; ch++){    

    fitmin = 3300;
    fitmax = 4000;    

    gh_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2+1]);
    gh_UnConnectedCh[ch]->Fit("pol1","","",1300,3300);
    sprintf(pltTit,"HG_%d",ch);
    Plot.TGraphPlotSetting(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2+1]);
    gl_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.TGraphPlotSetting(*gl_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2+1]);
    gtot_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.TGraphPlotSetting(*gtot_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    TF1* Linear_fit_h_Uncnct = gh_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l_Uncnct = gl_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot_Uncnct = gtot_UnConnectedCh[ch]->GetFunction("pol1");
    slope_h_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_l_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_tot_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);

    if(ch%50==0){
      sprintf(pltTit,"NcnctCh%d_InjCh%d_fit%d-%d",ch*2+1,Inj_ch,fitmin,fitmax);
      Plot.TGraphPlotSetting(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
			MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    }
  }

  //...Slope vs Uncnct Channel
  
  gh_Uncnct_Slope = new TGraph(NformatCH,UncnctID,slope_h_Uncnct);
  sprintf(pltTit,"UnconnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
  Plot.TGraphPlotSetting(*gh_Uncnct_Slope, pltTit, Xtit = "CH", Ytit = "P1",
		    MkSty = 22, MkClr = 2, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

  //...ConnectedCh

  for(int ch = 0; ch < NformatCH; ch++){
    fitmin = 0;
    fitmax = 1000;
    
    gh_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2]);
    gh_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"HG_%d",ch);
    Plot.TGraphPlotSetting(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2]);
    gl_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.TGraphPlotSetting(*gl_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2]);
    gtot_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.TGraphPlotSetting(*gtot_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    TF1* Linear_fit_h = gh_ConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l = gl_ConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot = gtot_ConnectedCh[ch]->GetFunction("pol1");
    slope_h[ch] = Linear_fit_h->GetParameter(1);
    slope_l[ch] = Linear_fit_l->GetParameter(1);
    slope_tot[ch] = Linear_fit_tot->GetParameter(1);

    CnctID[ch] = ch*2;

    if(ch%50==0){
      sprintf(pltTit,"NcnctCh%d_InjCh%d_fit%d-%d",ch*2+1,Inj_ch,fitmin,fitmax);
      Plot.TGraphPlotSetting(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
			MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    }
  }

  //...Slope vs cnct Channel
  
  gh_Cnct_Slope = new TGraph(NformatCH,CnctID,slope_h);
  sprintf(pltTit,"ConnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
  Plot.TGraphPlotSetting(*gh_Cnct_Slope, pltTit, Xtit = "Ch", Ytit = "P1",
		    MkSty = 22, MkClr = 3, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
  

  

  // ------------------------------ Plot HG_LG_TOT for Inj_ch all chip ------------------------------ //

  gh->Fit("pol1","","",0,200);
  TF1* Linear_fit_h_InjCh = gh->GetFunction("pol1");
  slope_h_InjCh = Linear_fit_h_InjCh->GetParameter(1);
  sprintf(pltTit,"Inj_CH_AllChip_Channel%dTS%d_HG",Inj_ch,MaxTS);
  Plot.TGraphPlotSetting(*gh, pltTit, Xtit = "DAC", Ytit = "ADC",
		    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
  
  gl->Fit("pol1","","",0,500);
  TF1* Linear_fit_l_InjCh = gl->GetFunction("pol1");
  slope_l_InjCh = Linear_fit_l_InjCh->GetParameter(1);
  sprintf(pltTit,"Inj_CH_AllChip_Channel%dTS%d_LG",Inj_ch,MaxTS);
  Plot.TGraphPlotSetting(*gl, pltTit, Xtit = "DAC", Ytit = "ADC",
		    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);  

  //------------------------------ Plotting HG_LG_TOT for Inj_ch per chip ------------------------------ //
  
  TLegend *legend = new TLegend(0.7,0.4,0.85,0.6);
  for(int ichip = 0; ichip < NCHIP; ichip++){
    gh_chip[ichip] = new TGraph(Nevents,dac_ctrl,ADC_H_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d",ichip);
    Plot.TGraphPlotSetting(*gh_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 24, MkClr = 1+ichip, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    
    gl_chip[ichip] = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d",ichip);
    Plot.TGraphPlotSetting(*gl_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 24, MkClr = 1+ichip, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_chip[ichip] = new TGraph(Nevents,dac_ctrl,TOT_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d",ichip);
    Plot.TGraphPlotSetting(*gtot_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
		      MkSty = 24, MkClr = 1+ichip, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    multig_InjCh_Chip_hltot->Add(gh_chip[ichip]);
    multig_InjCh_Chip_hltot->Add(gl_chip[ichip]);
    multig_InjCh_Chip_hltot->Add(gtot_chip[ichip]);
    legend->AddEntry(gh_chip[ichip],pltTit,"L");    
  }

  sprintf(pltTit,"HG_LG_TOT_InjCh%d",Inj_ch);
  Plot.TMultiGraphPlotSetting(*multig_InjCh_Chip_hltot, *legend, pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);


  
  //-------------------- Plotting the First Ring around the Inj_Ch -------------------- //
  
  TLegend *legendl = new TLegend(0.1,0.7,0.3,0.9);
  TLegend *legendh = new TLegend(0.1,0.7,0.3,0.9);
  for(int ichip = 0; ichip < NCHIP; ichip++){
    for(int cross_n = 0; cross_n < cross_num; cross_n++){
      if(cross_type_chip[cross_n][ichip] == true){
	gcross_h[cross_n] = new TGraph(Nevents,dac_ctrl,ADC_H_Cross_Chip[ichip][cross_n]);
	sprintf(pltTit,"CH %d High Gain",cross_ch_chip[cross_n][ichip]);
	Plot.TGraphPlotSetting(*gcross_h[cross_n], pltTit, Xtit = "DAC", Ytit = "ADC",
			  MkSty = 26, MkClr = cross_n+1, MkSize = 0.4, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
	legendh->AddEntry(gcross_h[cross_n],pltTit,"L");
	multig_cross_h->Add(gcross_h[cross_n]);

	gcross_l[cross_n] = new TGraph(Nevents,dac_ctrl,ADC_L_Cross_Chip[ichip][cross_n]);    
	sprintf(pltTit,"CH %d Low Gain",cross_ch_chip[cross_n][ichip]);
	Plot.TGraphPlotSetting(*gcross_l[cross_n], pltTit, Xtit = "DAC", Ytit = "ADC",
			  MkSty = 26, MkClr = cross_n+1, MkSize = 0.4, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
	legendl->AddEntry(gcross_l[cross_n],pltTit,"L");    
	multig_cross_l->Add(gcross_l[cross_n]);
      }
    }
    sprintf(pltTit,"FirstRing_around_Chip%dChannel%dLG",ichip,Inj_ch);
    Plot.TMultiGraphPlotSetting(*multig_cross_l, *legendl, pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
  }
  

  // ------------------------------ Cross Talk Ratio ------------------------------ //

  for(int ichip = 0; ichip < NCHIP; ichip++){
    gratioRing1_Injch_l[ichip] = new TGraph(Nevents,ADC_L_InjCh_Chip_double[ichip],Ratio_L_FirstRing_InjCh[ichip]);
    sprintf(pltTit,"Ratio_Ring1vsInjch%d_chip%d_L",Inj_ch,ichip);
    gratioRing1_Injch_l[ichip]->GetYaxis()->SetRangeUser(-0.1,0.1);
    Plot.TGraphPlotSetting(*gratioRing1_Injch_l[ichip], pltTit, Xtit = "LG_ADC", Ytit = "EFirstRing/ETotal",
			  MkSty = 26, MkClr = 1, MkSize = 0.4, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
  }
  

  //************************************************** Cross talk TH2Poly Plots **************************************************//

  int NNoisyCh = 7;
  int NoisyChannel[7] = {248,186,214,120,126,42,254};
  
  TH2Poly *polyh = new TH2Poly;
  TH2Poly *polyl = new TH2Poly;
  TH2Poly *polyInj = new TH2Poly;
  InitTH2Poly(*polyh);
  InitTH2Poly(*polyl);
  InitTH2Poly(*polyInj);
  
  for(int ch = 0; ch < NformatCH; ch++){
    float X, Y;
    bool NoisyBool = false;
    X = CHmap[ch].first;
    Y = CHmap[ch].second;
    if(Inj_ch%2!=1 && ch%32==Inj_ch/2){
      polyh->Fill(X,Y,0);
      polyl->Fill(X,Y,0);
      polyInj->Fill(X,Y,1);
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
  sprintf(pltTit,"Slope_HG_TS%dvsInjdac,Inj_ch=%d",MaxTS,Inj_ch);
  polyh->SetMaximum(0.2);
  Plot.TH2PolyPlotSetting(*polyh, pltTit, Xtit = "X[cm]", Ytit = "Y[cm]", Opt = "colztext", Stat = 0, Wait = 0, SavePlot = 0); 

  sprintf(pltTit,"Slope_LG_TS%dvsInjdac,Inj_ch=%d",MaxTS,Inj_ch);
  polyl->SetMaximum(0.05);
  Plot.TH2PolyPlotSetting(*polyl, pltTit, Xtit = "X[cm]", Ytit = "Y[cm]", Opt = "colztext", Stat = 0, Wait = 0, SavePlot = 0); 



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
  ifstream yamlFile(yamlFileName.c_str());
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

    /*    for(int header = 0; header < 13; header++) {getline(yamlFile,searchstr);}
    start = searchstr.find("'");
    searchstr = searchstr.substr(start+1,start+2);
    end = searchstr.find("',");
    searchstr = searchstr.erase(end);
    ModuleNumber = atoi(searchstr.c_str());
    cout << ModuleNumber << endl;
    */
    
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


void makePlots::Crosstalk(Int_t ichip, Int_t CH){

  TCanvas* c1 = new TCanvas();
  int cross_num = 6;
  int formatInj_Ch = CH/2+ichip*32;
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
    bool good_channel = true;
    while(abs(CHmap[ch].first-cross_posx[i]) > 1e-4 || abs(CHmap[ch].second-cross_posy[i]) > 1e-4){
      if(ch>127) {
	good_channel = false;
	break;
      }
      ch++;	    
      //cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
    }
    cross_ch[i] = (ch-ichip*32)*2;
    cross_type[i] = good_channel;
    cout << cross_ch[i] << endl;
    cout << cross_type[i] << endl;
  }
  getchar();
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

