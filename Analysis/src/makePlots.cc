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
  Crosstalk(Inj_ch);
  Plot.root_logon();
  if(Is_TB){
    sprintf(Plot.plotfolder_path,"plots/TBHexaboard/Injch_%d",Inj_ch);
  }
  else {
    sprintf(Plot.plotfolder_path,"plots/NTU_BarePCB/Injch_%d",Inj_ch);
  }  
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
  GainFactorReader();
  //outfile = new TFile("output.root","RECREATE");
  app = new TApplication("app",0,0);
}


//****************************************************************************************************//
//                                        PlotProducer                                                //
//****************************************************************************************************//



void makePlots::PlotProducer(){

  
  //==================== Define Parameters ====================

  char title[200];
  int TotalEntries = Chain1->GetEntries();
  int Nevents = TotalEntries/NCHIP;
  int AverageEvents = 0;
  int cross_num = 6;
  int MaxTS = 2; //choose this time sample to be the peak  
  
  int ADC_H_InjCh[Nevents], ADC_L_InjCh[Nevents], TOTS_InjCh[Nevents];
  int ADC_H_InjCh_Chip[NCHIP][Nevents], ADC_L_InjCh_Chip[NCHIP][Nevents], TOT_InjCh_Chip[NCHIP][Nevents];
  int ADC_H_FirstRing[NCHIP][Nevents], ADC_L_FirstRing[NCHIP][Nevents], TOT_FirstRing[NCHIP][Nevents];
  int ADC_H_ConnectedCh[NformatCH][Nevents],ADC_L_ConnectedCh[NformatCH][Nevents], TOT_ConnectedCh[NformatCH][Nevents];
  int ADC_H_AllCh[NCHANNEL][Nevents], ADC_L_AllCh[NCHANNEL][Nevents], TOT_AllCh[NCHANNEL][Nevents];
  int ADC_H_NoisyChannel[Nevents];
  double InjCh_MIP[NCHIP][Nevents], FirstRing_MIP[NCHIP][Nevents];
  double XTalkRatio_FirstRing[NCHIP][Nevents];
  double XTalkRatio_FirstRingChannels[NCHIP][cross_num][Nevents];
  double XTalkCoupling[NCHIP][NCH][Nevents], XTalkCoupling_Average[NCHIP][NCH];
  double ADC_L_InjCh_Chip_double[NCHIP][Nevents], Ratio_L_FirstRing_InjCh[NCHIP][event], ADC_H_InjCh_Chip_double[NCHIP][Nevents], Ratio_H_FirstRing_InjCh[NCHIP][event];
  double MIP[NCHIP][NCH][Nevents];
  int dac_ctrl[Nevents];
  double dac_ctrl_double[Nevents];
  bool cross_type_chip[cross_num][NCHIP];

  // Define Histograms
  TH1D *h_hgPedestal[NSCA];
  TH1D *h_lgPedestal[NSCA];

  for(int sca = 0; sca < NSCA; ++sca){
    char h_title[50];
    sprintf(h_title,"h_hgPedestal%d",sca);
    h_hgPedestal[sca] = new TH1D(h_title,"",100,0,500);
    sprintf(h_title,"h_lgPedestal%d",sca);
    h_lgPedestal[sca] = new TH1D(h_title,"",100,0,500);
  }
  /*
  TH1D *h = new TH1D("h","",100,150,250); //("title","",slice,star,end)
  TH1D *h_TOTS = new TH1D("h_TOTS","",100,5,500);
  TH1D *h_TOTF = new TH1D("h_TOTF","",100,1000,3000);
  TH1D *h_TOAR = new TH1D("h_TOAR","",100,1000,3000);
  TH1D *h_TOAF = new TH1D("h_TOAF","",100,1000,3000);
  */
  // Set Output Root File
  
  string Input_filename;
  Input_filename = input_RUN;
  int start = Input_filename.find_last_of("/");
  int end   = Input_filename.find(".root");
  string outf = Input_filename.substr(start+1,end-start-1);

  sprintf(title,"output_root/%s.root",outf.c_str());
  TFile *outfile = new TFile(title,"recreate");
  //TTree *outT1 = new TTree("Rechit_var","Rechit_var");

  sprintf(title,"hgPedestal.txt");


  //==================== Initialize ====================
  for(int ichip = 0; ichip < NCHIP; ichip++){
    for(int ich = 0; ich < NCH; ich++){
      XTalkCoupling_Average[ichip][ich] = 0;
    }
  }

  
  for (int ev = 0; ev < Nevents; ev++){
    ADC_H_InjCh[ev] = 0;
    ADC_L_InjCh[ev] = 0;
    dac_ctrl[ev] = ev;
    for(int ichip = 0; ichip < NCHIP; ichip++){
      ADC_H_FirstRing[ichip][ev] = 0;
      ADC_L_FirstRing[ichip][ev] = 0;
      TOT_FirstRing[ichip][ev] = 0;
      for (int icross = 0; icross < cross_num; icross++){
	XTalkRatio_FirstRingChannels[ichip][icross][ev] = 0;
      }
    }
  }
    
  
  //==================== Loop over the events ====================
   
  for(int entry = 0; entry < TotalEntries ; ++entry){
    
    if(entry%1000==0){ cout << "Now Processing entry = " << entry << endl; }    
    Chain1 -> GetEntry(entry);
    dac_ctrl[event] = dacinj;
    dac_ctrl_double[event] = dacinj;    
    
    
    int TS0_sca, MaxTS_sca;
    for(int sca = 0 ; sca < NSCA ; sca++) {
      TS[sca] = timesamp[sca];
      if (timesamp[sca] == 0) { TS0_sca = sca ; }
      if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
    }

    // -------------------- Pedestal Analysis -------------------- //
    for(int sca = 0; sca < NSCA; ++sca){
	  if ( chip != 3 ) continue;
      h_hgPedestal[sca]->Fill(hg[sca][52]); // Fill histogram with TS0 readout 
      h_lgPedestal[sca]->Fill(lg[sca][52]);
    }

    
    // -------------------- Injection & Cross Talk Analysis -------------------- //

    // Filling hg, lg data (with 13 SCA)	
    ADC_H_InjCh[event] = hg[MaxTS_sca][Inj_ch]; // Filling the Inj_ch
    ADC_L_InjCh[event] = lg[MaxTS_sca][Inj_ch];

    ADC_H_InjCh_Chip[chip][event] = hg[MaxTS_sca][Inj_ch];
    ADC_L_InjCh_Chip[chip][event] = lg[MaxTS_sca][Inj_ch];
    ADC_H_InjCh_Chip_double[chip][event] = ( hg[MaxTS_sca][Inj_ch]  );
    ADC_L_InjCh_Chip_double[chip][event] = ( lg[MaxTS_sca][Inj_ch]  );
    TOT_InjCh_Chip[chip][event] = tot_slow[Inj_ch];

    
    
    // Calculate energy converted into mip 
    for(int ich = 0; ich < NCH; ich++){
      int hg_sig = hg[MaxTS_sca][ich] - hg[TS0_sca][ich];
      int lg_sig = lg[MaxTS_sca][ich] - lg[TS0_sca][ich];
      int tot_sig = tot_slow[ich];

      
      if( hg_sig < HGTP[chip][ich]){
	MIP[chip][ich][event] = hg_sig * ADC2MIP;
      }
      else{
	//if( lg_sig < LGTP[chip][ich]){
	if( lg_sig < LGTP_default){
	  MIP[chip][ich][event] = ( lg_sig * LG2HG_Conversion[chip][ich] * ADC2MIP);
	}
	else{
	  MIP[chip][ich][event] = ( (tot_sig - TOTOffSet[chip][ich]) * TOT2LG_Conversion[chip][ich] * LG2HG_Conversion[chip][ich] * ADC2MIP); 
	}
      }
    }

    
    for(int ch = 0; ch < 32; ch++){
      ADC_H_ConnectedCh[ch+chip*32][event] = hg[MaxTS_sca][ch*2]; // Filling all the connected channels
      ADC_L_ConnectedCh[ch+chip*32][event] = lg[MaxTS_sca][ch*2];
      ADC_H_AllCh[ch+chip*64][event] = hg[MaxTS_sca][ch];
      ADC_H_AllCh[ch+chip*64][event] = hg[MaxTS_sca][ch];
    }

    // XTalk coupling calculation
    for(int ich = 0; ich < NCH; ich++){
      XTalkCoupling[chip][ich][event] = MIP[chip][ich][event] / MIP[chip][Inj_ch][event];
      if( event>50 && event<=700 ){
	XTalkCoupling_Average[chip][ich] += XTalkCoupling[chip][ich][event];
	AverageEvents++;
      }
    }


    // FirstRing Analsis
    for(int icross = 0; icross < cross_num; icross++){
      if(cross_type[chip][icross]){ 
	int xchip = cross_ch_FirstRing[chip][icross] % NCH;
	int xch = cross_ch_FirstRing[chip][icross]- (xchip * NCH);
	double xMIP = MIP[xchip][xch][event], injMIP = MIP[chip][Inj_ch][event];

	FirstRing_MIP[chip][event] += xMIP;
	XTalkRatio_FirstRingChannels[chip][icross][event] = xMIP / injMIP;
      }
    }
    XTalkRatio_FirstRing[chip][event] = FirstRing_MIP[chip][event]/(MIP[chip][Inj_ch][event] + FirstRing_MIP[chip][event]);

	

    // Filling tot, toa data (without SCA)
    for(int ch = 0; ch < NformatCH/NCHIP; ch++){
      TOT_ConnectedCh[ch+chip*32][event] = tot_slow[ch*2];
    }
  }

  //... ==================== End of Loop ==================== ...

  for(int ichip = 0; ichip < NCHIP; ichip++){
    for(int ich = 0; ich < NCH; ich++){
      XTalkCoupling_Average[ichip][ich] /= (AverageEvents/NCHANNEL);
    }
  }
  
  //... ==================== Fit & Plots ==================== ...
  

  // Define Plotting Paramter
  
  char pltTit[200], leg[50], img_title[50];
  string Xtit, Ytit, Opt;
  int MkSty, MkClr, LClr, fitmin, fitmax;
  float MkSize, LWid;
  bool Stat, Wait, SavePlot;


  // Define Fitting Parameter
  /*
  double slope_h[NformatCH], slope_l[NformatCH], slope_tot[NformatCH];
  double slope_h_Uncnct[NformatCH], slope_l_Uncnct[NformatCH], slope_tot_Uncnct[NformatCH];
  double slope_h_InjCh, slope_l_InjCh, slope_h_chip[NCHIP], slope_l_chip[NCHIP];
  double CnctID[NformatCH], UncnctID[NformatCH];

  */


  // Plots!!!!!

  TMultiGraph* multig_XRatioRing1 = new TMultiGraph();
  TLegend* legend_Ring1 = new TLegend(0.5, 0.92, 0.9, 0.98);
  legend_Ring1->SetNColumns(4);
  TGraph* gXTalkRatioRing1 = new TGraph();

  for(int ichip = 0; ichip < NCHIP; ichip++){
    
    // ------------------------------ MIP Conversion Plot ------------------------------ //
    TGraph* gInjch_mip = new TGraph(Nevents,dac_ctrl_double,MIP[ichip][Inj_ch]);
    sprintf(pltTit,"InjCh%d_chip%d_MIP",Inj_ch,ichip);
    Plot.GStd(*gInjch_mip, pltTit, Xtit = "DAC", Ytit = "MIP", Opt = "AP", Wait = 0, SavePlot = 0);
    
    // ------------------------------ Cross Talk Ratio ------------------------------ //
    gXTalkRatioRing1 = new TGraph(Nevents,MIP[ichip][Inj_ch],XTalkRatio_FirstRing[ichip]);
    sprintf(pltTit,"Chip%d", ichip);
    Plot.MultiAdd(*multig_XRatioRing1, *gXTalkRatioRing1, *legend_Ring1, pltTit, MkSty = 25, MkClr = ichip, MkSize = 0.5);

    // -----------------------------  Cross Talk Ratio based on each FirstRing channel ------------------------------ //
    
    TMultiGraph* multig_XRatioRing1Channels = new TMultiGraph();
    TLegend* legend_Ring1Channels = new TLegend(0.53, 0.75, 0.88, 0.88);
    legend_Ring1Channels->SetNColumns(2);
    
    for(int icross = 0; icross < cross_num; icross++){
      if( cross_type[ichip][icross] == true) {
	TGraph* gXTalkRatioRing1Channels = new TGraph(Nevents, MIP[ichip][Inj_ch], XTalkRatio_FirstRingChannels[ichip][icross]);
	sprintf(pltTit,"Pos_%d, Ch_%d", icross, cross_ch_FirstRing[ichip][icross]);
	Plot.MultiAdd(*multig_XRatioRing1Channels, *gXTalkRatioRing1Channels, *legend_Ring1Channels, pltTit, MkSty = 25, MkClr = icross, MkSize = 0.5);
      }
    }
    multig_XRatioRing1Channels->SetMaximum(0.05);
    multig_XRatioRing1Channels->SetMinimum(0.);
    sprintf(pltTit,"XTalkCoupling_FirstRing_InjCh%d_Chip%d",Inj_ch, ichip);
    Plot.Multi(*multig_XRatioRing1Channels, *legend_Ring1Channels, pltTit, Xtit = "EInj[MIP]", Ytit = "EFirstRing / EInj", Opt = "AP", Wait = 0, SavePlot = 0);
  }

  multig_XRatioRing1->SetMaximum(0.2);
  multig_XRatioRing1->SetMinimum(-0.1);
  sprintf(pltTit,"XTalkCoupling_FirstRing_InjCh%d",Inj_ch);
  Plot.Multi(*multig_XRatioRing1, *legend_Ring1, pltTit, Xtit = "EInj[MIP]", Ytit = "EFirstRing / EInj", Opt = "AP", Wait = 0, SavePlot = 0);

  TGraph* gexample = new TGraph(Nevents,MIP[0][Inj_ch],XTalkRatio_FirstRingChannels[0][0]);
  sprintf(pltTit,"Chip%d_Ch%d_InjCh%d",0,0,Inj_ch);
  gexample->GetYaxis()->SetRangeUser(-0.,0.05);
  Plot.GStd(*gexample, pltTit, Xtit = "EInj[MIP]", Ytit = "E / EInj", Opt = "AP", Wait = 0, SavePlot = 0);
    
  
  // 2D Graph of XTalk coupling

  int NNoisy = 8;
  int NoisyChannel[8] = {248,186,214,120,126,42,254,190};

  TH2Poly *poly = new TH2Poly;
  InitTH2Poly(*poly);
  poly->SetMinimum(-0.005);
  for(int ichip = 0; ichip < NCHIP; ichip++){
    for(int ich = 0; ich < NCH; ich++){
      float X, Y;
      int forCH = ichip*(NCH/2) + ich/2;
      bool NoisyBool = false;
      X = CHmap[forCH].first;
      Y = CHmap[forCH].second;
      if(ich == Inj_ch){
	poly->Fill(X,Y,-2);
      }
      if( ich != Inj_ch && ich%2 == 0){
	/*	for(int iNoisy = 0; iNoisy < NNoisy; iNoisy++){
	  if(forCH == NoisyChannel[iNoisy]/2) {NoisyBool = true;}
	  }*/
	if(!NoisyBool){
	  //poly->Fill(X,Y,XTalkCoupling_Average[ichip][ich]);
	  poly->Fill(X,Y,forCH);
	}
      }
    }
  }
  sprintf(pltTit,"InjCh%d",Inj_ch);
  Plot.Poly(*poly, pltTit, Xtit = "X[cm]", Ytit = "Y[cm]", Opt = "colz", Stat = 0, Wait = 0, SavePlot = 0);

  // 1D Graph  


  // Unconnect channels
  
  double LabelUncnct[NformatCH];
  double XTalkCoupling_Uncnct[NformatCH];
  for(int ichip = 0; ichip < NCHIP; ichip++){
    for(int ich = 1; ich < NCH; ich = ich+2){
      LabelUncnct[(ich+ichip*NCH)/2] = ich + ichip*NCH;
      XTalkCoupling_Uncnct[(ich+ichip*NCH)/2] = XTalkCoupling_Average[ichip][ich];
    }
  }
  TGraph* gXTalkCoupling_Uncnct = new TGraph(NformatCH, LabelUncnct , XTalkCoupling_Uncnct);
  sprintf(pltTit,"XTalkCoupling_Uncnct_InjCh%d",Inj_ch);
  Plot.GStd(*gXTalkCoupling_Uncnct, pltTit, Xtit = "UnconnectCh ID", Ytit = "E / EInj", Opt = "AP", Wait = 0, SavePlot = 0);

  TGraph* gXTalkCoupling_Uncnct39 = new TGraph(Nevents, MIP[0][Inj_ch], MIP[0][39]);
  sprintf(pltTit,"Uncnct39");
  Plot.GStd(*gXTalkCoupling_Uncnct39, pltTit, Xtit = "EInj[MIP]", Ytit = "E / EInj", Opt = "AP", Wait = 0, SavePlot = 0);



  //------------------------------ Plotting HG_LG_TOT for Inj_ch per chip ------------------------------ //

  TGraph **gh_chip = new TGraph*[NCHIP];
  TGraph **gl_chip = new TGraph*[NCHIP];
  TGraph **gtot_chip = new TGraph*[NCHIP];
  
  for(int ichip = 0; ichip < NCHIP; ichip++){
    TLegend *legend = new TLegend(0.7,0.4,0.85,0.6);
    TMultiGraph *multig_InjCh_Chip_hltot = new TMultiGraph();
    gh_chip[ichip] = new TGraph(Nevents,dac_ctrl,ADC_H_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d_inj%d",ichip,Inj_ch);
    Plot.G(*gh_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
	   MkSty = 24, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    legend->AddEntry(gh_chip[ichip],"HG","L");
    
    gl_chip[ichip] = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d_inj%d",ichip,Inj_ch);
    Plot.G(*gl_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
	   MkSty = 24, MkClr = 2, MkSize = 0.5, LClr = 2, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    legend->AddEntry(gl_chip[ichip],"LG","L");

    gtot_chip[ichip] = new TGraph(Nevents,dac_ctrl,TOT_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d_inj%d",ichip,Inj_ch);
    Plot.G(*gtot_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
	   MkSty = 24, MkClr = 3, MkSize = 0.5, LClr = 3, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    legend->AddEntry(gtot_chip[ichip],"TOT","L");

    multig_InjCh_Chip_hltot->Add(gh_chip[ichip]);
    multig_InjCh_Chip_hltot->Add(gl_chip[ichip]);
    multig_InjCh_Chip_hltot->Add(gtot_chip[ichip]);
    sprintf(pltTit,"HG_LG_TOT_InjCh%d_Chip%d",Inj_ch,ichip);
    Plot.Multi(*multig_InjCh_Chip_hltot, *legend, pltTit, Xtit = "Injection DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);
    delete multig_InjCh_Chip_hltot;
    delete legend;
  }

  // ------------------------------ Pedestal Plots ------------------------------ //

  for(int sca = 0; sca < NSCA; ++sca){
    sprintf(pltTit,"hgPedestal[%d]",sca);
    Plot.HStd(*h_hgPedestal[sca], pltTit, Xtit = "ADC", Ytit = "", Wait = 0, SavePlot = 0);
	h_hgPedestal[sca]->Write();

	sprintf(pltTit,"lgPedestal[%d]",sca);
    Plot.HStd(*h_lgPedestal[sca], pltTit, Xtit = "ADC", Ytit = "", Wait = 0, SavePlot = 0);
	h_lgPedestal[sca]->Write();

  }


  outfile->Write();
  outfile->Close();


  // ------------------------------ Fit ------------------------------ //

  // UnconnectedCh
  /*
    for(int ch = 0; ch < NformatCH; ch++){    

    fitmin = 1300;
    fitmax = 3300;    

    gh_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2+1]);
    gh_UnConnectedCh[ch]->Fit("pol1","","",1300,3300);
    sprintf(pltTit,"HG_%d",ch);
    Plot.G(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2+1]);
    gl_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.G(*gl_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2+1]);
    gtot_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.G(*gtot_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    TF1* Linear_fit_h_Uncnct = gh_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l_Uncnct = gl_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot_Uncnct = gtot_UnConnectedCh[ch]->GetFunction("pol1");
    slope_h_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_l_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_tot_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);

    if(ch%50==0){
    sprintf(pltTit,"NcnctCh%d_InjCh%d_fit%d-%d",ch*2+1,Inj_ch,fitmin,fitmax);
    Plot.G(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);
    }
    }
  
    //...Slope vs Uncnct Channel
  
    gh_Uncnct_Slope = new TGraph(NformatCH,UncnctID,slope_h_Uncnct);
    sprintf(pltTit,"UnconnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
    Plot.G(*gh_Uncnct_Slope, pltTit, Xtit = "CH", Ytit = "P1",
    MkSty = 22, MkClr = 2, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);


    //... ConnectedCh

    for(int ch = 0; ch < NformatCH; ch++){
    fitmin = 1000;
    fitmax = 4000;
    
    gh_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2]);
    gh_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"HG_%d",ch);
    Plot.G(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2]);
    gl_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.G(*gl_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2]);
    gtot_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.G(*gtot_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
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
    Plot.G(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);
    }
    }

    //...Slope vs cnct Channel
  
    gh_Cnct_Slope = new TGraph(NformatCH,CnctID,slope_h);
    sprintf(pltTit,"ConnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
    Plot.G(*gh_Cnct_Slope, pltTit, Xtit = "Ch", Ytit = "P1",
    MkSty = 22, MkClr = 3, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);


    //------------------------------ Change fitting Range ------------------------------//

    //... UnconncetedCh

    for(int ch = 0; ch < NformatCH; ch++){    

    fitmin = 3300;
    fitmax = 4000;    

    gh_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2+1]);
    gh_UnConnectedCh[ch]->Fit("pol1","","",1300,3300);
    sprintf(pltTit,"HG_%d",ch);
    Plot.G(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2+1]);
    gl_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.G(*gl_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_UnConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2+1]);
    gtot_UnConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.G(*gtot_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    TF1* Linear_fit_h_Uncnct = gh_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_l_Uncnct = gl_UnConnectedCh[ch]->GetFunction("pol1");
    TF1* Linear_fit_tot_Uncnct = gtot_UnConnectedCh[ch]->GetFunction("pol1");
    slope_h_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_l_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);
    slope_tot_Uncnct[ch] = Linear_fit_h_Uncnct->GetParameter(1);

    if(ch%50==0){
    sprintf(pltTit,"NcnctCh%d_InjCh%d_fit%d-%d",ch*2+1,Inj_ch,fitmin,fitmax);
    Plot.G(*gh_UnConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);
    }
    }

    //...Slope vs Uncnct Channel
  
    gh_Uncnct_Slope = new TGraph(NformatCH,UncnctID,slope_h_Uncnct);
    sprintf(pltTit,"UnconnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
    Plot.G(*gh_Uncnct_Slope, pltTit, Xtit = "CH", Ytit = "P1",
    MkSty = 22, MkClr = 2, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);

    //...ConnectedCh

    for(int ch = 0; ch < NformatCH; ch++){
    fitmin = 0;
    fitmax = 1000;
    
    gh_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_H_AllCh[ch*2]);
    gh_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"HG_%d",ch);
    Plot.G(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gl_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,ADC_L_AllCh[ch*2]);
    gl_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"LG_%d",ch);
    Plot.G(*gl_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_ConnectedCh[ch] = new TGraph(Nevents,dac_ctrl,TOT_AllCh[ch*2]);
    gtot_ConnectedCh[ch]->Fit("pol1","","",fitmin,fitmax);
    sprintf(pltTit,"TOT_%d",ch);
    Plot.G(*gtot_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
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
    Plot.G(*gh_ConnectedCh[ch], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);
    }
    }

    //...Slope vs cnct Channel
  
    gh_Cnct_Slope = new TGraph(NformatCH,CnctID,slope_h);
    sprintf(pltTit,"ConnectedCh_InjCh%d_fit%d-%d",Inj_ch,fitmin,fitmax);
    Plot.G(*gh_Cnct_Slope, pltTit, Xtit = "Ch", Ytit = "P1",
    MkSty = 22, MkClr = 3, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);
  

  

    // ------------------------------ Plot HG_LG_TOT for Inj_ch all chip ------------------------------ //

    gh->Fit("pol1","","",0,200);
    TF1* Linear_fit_h_InjCh = gh->GetFunction("pol1");
    slope_h_InjCh = Linear_fit_h_InjCh->GetParameter(1);
    sprintf(pltTit,"Inj_CH_AllChip_Channel%dTS%d_HG",Inj_ch,MaxTS);
    Plot.G(*gh, pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
  
    gl->Fit("pol1","","",0,500);
    TF1* Linear_fit_l_InjCh = gl->GetFunction("pol1");
    slope_l_InjCh = Linear_fit_l_InjCh->GetParameter(1);
    sprintf(pltTit,"Inj_CH_AllChip_Channel%dTS%d_LG",Inj_ch,MaxTS);
    Plot.G(*gl, pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 7, MkClr = 1, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);  

    //------------------------------ Plotting HG_LG_TOT for Inj_ch per chip ------------------------------ //
  
    TLegend *legend = new TLegend(0.7,0.4,0.85,0.6);
    for(int ichip = 0; ichip < NCHIP; ichip++){
    gh_chip[ichip] = new TGraph(Nevents,dac_ctrl,ADC_H_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d",ichip);
    Plot.G(*gh_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 24, MkClr = 1+ichip, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    
    gl_chip[ichip] = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d",ichip);
    Plot.G(*gl_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 24, MkClr = 1+ichip, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    gtot_chip[ichip] = new TGraph(Nevents,dac_ctrl,TOT_InjCh_Chip[ichip]);
    sprintf(pltTit,"Chip%d",ichip);
    Plot.G(*gtot_chip[ichip], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 24, MkClr = 1+ichip, MkSize = 0.5, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);

    multig_InjCh_Chip_hltot->Add(gh_chip[ichip]);
    multig_InjCh_Chip_hltot->Add(gl_chip[ichip]);
    multig_InjCh_Chip_hltot->Add(gtot_chip[ichip]);
    legend->AddEntry(gh_chip[ichip],pltTit,"L");    
    }

    sprintf(pltTit,"HG_LG_TOT_InjCh%d",Inj_ch);
    Plot.Multi(*multig_InjCh_Chip_hltot, *legend, pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);


  
    //-------------------- Plotting the First Ring around the Inj_Ch -------------------- //
  
    TLegend *legendl = new TLegend(0.1,0.7,0.3,0.9);
    TLegend *legendh = new TLegend(0.1,0.7,0.3,0.9);
    for(int ichip = 0; ichip < NCHIP; ichip++){
    for(int cross_n = 0; cross_n < cross_num; cross_n++){
    if(cross_type_chip[cross_n][ichip] == true){
    gcross_h[cross_n] = new TGraph(Nevents,dac_ctrl,ADC_H_Cross_Chip[ichip][cross_n]);
    sprintf(pltTit,"CH %d High Gain",cross_ch_chip[cross_n][ichip]);
    Plot.G(*gcross_h[cross_n], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 26, MkClr = cross_n+1, MkSize = 0.4, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    legendh->AddEntry(gcross_h[cross_n],pltTit,"L");
    multig_cross_h->Add(gcross_h[cross_n]);

    gcross_l[cross_n] = new TGraph(Nevents,dac_ctrl,ADC_L_Cross_Chip[ichip][cross_n]);    
    sprintf(pltTit,"CH %d Low Gain",cross_ch_chip[cross_n][ichip]);
    Plot.G(*gcross_l[cross_n], pltTit, Xtit = "DAC", Ytit = "ADC",
    MkSty = 26, MkClr = cross_n+1, MkSize = 0.4, LClr = 1, LWid = 4, Opt = "AP", Stat = 1, Wait = 0, SavePlot = 0);
    legendl->AddEntry(gcross_l[cross_n],pltTit,"L");    
    multig_cross_l->Add(gcross_l[cross_n]);
    }
    }
    sprintf(pltTit,"FirstRing_around_Chip%dChannel%dLG",ichip,Inj_ch);
    Plot.Multi(*multig_cross_l, *legendl, pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Stat = 1, Wait = 0, SavePlot = 1);
    }



    //========================================//
    */
  
  
  //************************************************** Cross talk TH2Poly Plots **************************************************//
  /*  
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
  Plot.Poly(*polyh, pltTit, Xtit = "X[cm]", Ytit = "Y[cm]", Opt = "colztext", Stat = 0, Wait = 0, SavePlot = 0); 

  sprintf(pltTit,"Slope_LG_TS%dvsInjdac,Inj_ch=%d",MaxTS,Inj_ch);
  polyl->SetMaximum(0.05);
  Plot.Poly(*polyl, pltTit, Xtit = "X[cm]", Ytit = "Y[cm]", Opt = "colztext", Stat = 0, Wait = 0, SavePlot = 0); 
  
*/

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
  if(TS[sca] == MaxTS ){
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


  int start = input_RUN.find_last_of("/");
  int end = input_RUN.find("_pedestal.root");
  string f_substr = input_RUN.substr(start+1,end-start-1);
  string rootFileName(f_substr);
  end = input_RUN.find("ana_output");
  cout << end << endl;
  f_substr = input_RUN.substr(0,end-1);
  string yamlPath(f_substr);
  char yamlFileName[100];
  sprintf(yamlFileName,"%s/yaml/%s.yaml",yamlPath.c_str(), rootFileName.c_str());
  
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
    cout << "InjCh = " << Inj_ch << endl;

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

void makePlots::GainFactorReader(){

  string GainFileName;
  string TB_GainFactors("src_txtfile/TPro_fittingoutput.txt");
  if(Is_TB){
    GainFileName = TB_GainFactors;}
  else{
    string GainFileName("src_txtfile/TPro_NTU_Inj.txt");}
  
  ifstream GainFile(GainFileName.c_str());
  string line;
  char tmp[50];
  int ichip, ich;
  if(!GainFile.is_open()){
    cout << "Did not find GainFactor file " << GainFileName
	 << ".\n Take this run's GainFactor as default & calculate the Gainfactors " << endl;
    for(ichip = 0; ichip < NCHIP; ichip++){
      for(ich = 0; ich < NCH; ich++){
	HGTP[ichip][ich] = 1500;
	LG2HG_Conversion[ichip][ich] = 8.5;
	LGTP[ichip][ich] = 900;
	TOT2LG_Conversion[ichip][ich] = 3.8;
	TOTOffSet[ichip][ich] = 180;
      }
    }
    Gain_factor_producer();
  }
  
  if(GainFile.is_open()){
    getline(GainFile,line);
    while(!GainFile.eof()){
      GainFile >> tmp >> tmp >> ichip >> ich >> tmp;
      GainFile >> HGTP[ichip][ich] >> LG2HG_Conversion[ichip][ich] >> LGTP[ichip][ich] >> TOT2LG_Conversion[ichip][ich] >> TOTOffSet[ichip][ich];
      getline(GainFile,line);
    }
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
    if(ev % 50 != 0) continue;
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
  for(int ichip = 0; ichip < NCHIP; ichip++){
    
    int cross_num = 6;
    int formatInj_Ch = CH/2+ichip*32;
    float Xdist = 0.974452;
    float Ydist = 0.5626;
    float cross_posx[cross_num], cross_posy[cross_num];
    float X = CHmap[formatInj_Ch].first;
    float Y = CHmap[formatInj_Ch].second;
  
    cross_posx[0] = X - Xdist;
    cross_posy[0] = Y + Ydist;
    cross_posx[1] = X;
    cross_posy[1] = Y + 2*Ydist;
    cross_posx[2] = X + Xdist;
    cross_posy[2] = Y + Ydist;
    cross_posx[3] = X + Xdist;
    cross_posy[3] = Y - Ydist;
    cross_posx[4] = X;
    cross_posy[4] = Y - 2*Ydist;
    cross_posx[5] = X - Xdist;
    cross_posy[5] = Y - Ydist;
  
    cout << "Chip" << ichip << " FirstRing Channels = " ;
    for(int icross = 0; icross < 6; icross++){
      //cout << cross_posx[i] << " " << cross_posy[i] << endl;
      int ch = 0;
      bool good_channel = true;
      while(abs(CHmap[ch].first-cross_posx[icross]) > 1e-4 || abs(CHmap[ch].second-cross_posy[icross]) > 1e-4){
	if( ch > 127) {
	  good_channel = false;
	  break;
	}
	ch++;	    
	//cout << ch << " "<<abs(CHmap[ch].first - cross_posx[i])<< " " << abs(CHmap[ch].second - cross_posy[i])<< endl;
	/*
	forCH = chip*32 + cross_ch[i]/2;
	X = CHmap[forCH].first;
	Y = CHmap[forCH].second;	      
	poly->Fill(X,Y,cross_ch[i]);
	*/
      
      }
      cross_ch_FirstRing[ichip][icross] = ch * 2;
      cross_type[ichip][icross] = good_channel;      
      cout << cross_ch_FirstRing[ichip][icross] << " " ;
    }
    cout << endl;
  }
  
  TH2Poly *poly = new TH2Poly;
  InitTH2Poly(*poly);
  /*  
  for(int i=0; i<6; i++){
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
      if (timesamp[sca] == 0) { TS0_sca = sca ; }
      if (timesamp[sca] == MaxTS) { MaxTS_sca = sca ; }
    }
    
    
    ADC_H_InjCh_Chip[chip][event] = ( hg[MaxTS_sca][Inj_ch] - hg[TS0_sca][Inj_ch] );
    ADC_L_InjCh_Chip[chip][event] = ( lg[MaxTS_sca][Inj_ch] - lg[TS0_sca][Inj_ch] );
    TOT_InjCh_Chip[chip][event] = tot_slow[Inj_ch];


    if(ADC_H_InjCh_Chip[chip][event] > HGTP[chip][Inj_ch] && HGTP_flag[chip] == false){
      HGTP_flag[chip] = true;
      HGLGfitmax[chip] = ADC_L_InjCh_Chip[chip][event];
    }
    if(ADC_L_InjCh_Chip[chip][event] > LGTP[chip][Inj_ch] && LGTP_flag[chip] == false){
      LGTP_flag[chip] = true;
      TOTLGfitmax[chip] = TOT_InjCh_Chip[chip][event];
    }
  }

  
  //==================== End Loop ====================

  //...

  //==================== Plots ====================

  char pltTit[100];
  string Xtit, Ytit, Opt;
  int MkSty, MkClr, LClr, fitmin, fitmax;
  float MkSize, LWid;
  bool Stat, Wait, SavePlot;
  TCanvas *c1 = new TCanvas();

  TGraph** gh = new TGraph*[NCHIP];
  TGraph** gl = new TGraph*[NCHIP];
  TGraph** gtot = new TGraph*[NCHIP];
  TGraph** LG2HG = new TGraph*[NCHIP];
  TGraph** TOT2LG = new TGraph*[NCHIP];

  for(int ichip = 0; ichip < NCHIP; ichip++){
    gh[ichip] = new TGraph(Nevents,dac_ctrl,ADC_H_InjCh_Chip[ichip]);
    gl[ichip] = new TGraph(Nevents,dac_ctrl,ADC_L_InjCh_Chip[ichip]);
    gtot[ichip] = new TGraph(Nevents,dac_ctrl,TOT_InjCh_Chip[ichip]);
    LG2HG[ichip] = new TGraph(Nevents,ADC_L_InjCh_Chip[ichip],ADC_H_InjCh_Chip[ichip]);
    TOT2LG[ichip] = new TGraph(Nevents,TOT_InjCh_Chip[ichip],ADC_L_InjCh_Chip[ichip]);
    
    LG2HG[ichip]->Fit("pol1","","",fitmin = 0,fitmax = HGLGfitmax[ichip]);
    //TOT2LG[ichip]->Fit("pol1","","",fitmin = TOTOffSet,fitmax = TOTLGfitmax[ichip]);
    TOT2LG[ichip]->Fit("pol1","","",fitmin = 200,fitmax = 300);
    
    TF1* Linear_fit_LG2HG = LG2HG[ichip]->GetFunction("pol1");
    TF1* Linear_fit_TOT2LG = TOT2LG[ichip]->GetFunction("pol1");
    LG2HG_Conversion[ichip][Inj_ch] = Linear_fit_LG2HG->GetParameter(1);
    TOT2LG_Conversion[ichip][Inj_ch] = Linear_fit_TOT2LG->GetParameter(1);
    TOTOffSet[ichip][Inj_ch] = -Linear_fit_TOT2LG->GetParameter(0)/Linear_fit_TOT2LG->GetParameter(1);
    cout << LG2HG_Conversion[ichip][Inj_ch] << " " <<  TOT2LG_Conversion[ichip][Inj_ch] << " " << TOTOffSet[ichip][Inj_ch] << endl;
    /*
      TOT2LG[ichip]->Draw("AP");
      c1->Update();
      gPad->WaitPrimitive();
    */
    sprintf(pltTit,"HG_Chip%d",ichip);
    Plot.GStd(*gh[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);

    sprintf(pltTit,"LG_Chip%d",ichip);
    Plot.GStd(*gl[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);

    sprintf(pltTit,"TOT_Chip%d",ichip);
    Plot.GStd(*gtot[ichip], pltTit, Xtit = "DAC", Ytit = "ADC", Opt = "AP", Wait = 0, SavePlot = 1);
    
    sprintf(pltTit,"LG2HG_Chip%d",ichip);
    Plot.GStd(*LG2HG[ichip], pltTit, Xtit = "LG", Ytit = "HG", Opt = "AP", Wait = 0, SavePlot = 1);
    
    sprintf(pltTit,"TOT2LG_Chip%d",ichip);
    Plot.GStd(*TOT2LG[ichip], pltTit, Xtit = "TOT", Ytit = "LG", Opt = "AP", Wait = 0, SavePlot = 1);
  }
}



