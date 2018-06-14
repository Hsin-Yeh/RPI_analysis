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

ClassImp(hit)
ClassImp(hitcollection)
//Constructo
makePlots::makePlots(TChain* inchain):fChain(inchain)
{
  HITS = new hitcollection;
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
}

void makePlots::Loop(){

  Init();
  //P_and_N(0,1);
  //read_P_and_N("ped_result/Module_1_RUN_300318_0527");

  /*
  int BeforeHG[4][64][NSCA], BeforeLG[4][64][NSCA];
  int AfterHG[4][64][NSCA], AfterLG[4][64][NSCA];
  TH1D* hpedcheck = new TH1D("pedcheck","Low Gain Pedestal Check (End - Start)",20,-10,10);

  for(int chip=0; chip<4; ++chip){
    for(int ch=0; ch<64; ++ch){ 
      for(int sca = 0; sca < NSCA; ++sca){
	BeforeHG[chip][ch][sca] = avg_HG[chip][ch][sca]; 
	BeforeLG[chip][ch][sca] = avg_LG[chip][ch][sca];
      }
    }
  }

  read_P_and_N("ped_result/Module_1_RUN_300318_0659");

  for(int chip=0; chip<4; ++chip){
    for(int ch=0; ch<64; ++ch){ 
      for(int sca = 0; sca < NSCA; ++sca){
	avg_HG[chip][ch][sca] -= BeforeHG[chip][ch][sca]; 
	avg_LG[chip][ch][sca] -= BeforeLG[chip][ch][sca];
	AfterHG[chip][ch][sca] = avg_HG[chip][ch][sca]; 
	AfterLG[chip][ch][sca] = avg_LG[chip][ch][sca];
	hpedcheck->Fill(AfterLG[chip][ch][sca]);
      }
    }
  }
  */
  
  

  

  /*  for(int chip = 0; chip < NCHIP; chip++){
      for(int ch = 0; ch < 64; ch++){
      for(int sca = 0; sca < NSCA ; ++sca){
      cout << avg_HG[chip][ch][sca] << endl;}}}
  */
  //  avg_HG[chip][ch][sca]

  //==================== Call the Parameters ====================
  
  app = new TApplication("app",0,0);
  TCanvas *c1 = new TCanvas;
  //  TCanvas *cmulti = new TCanvas;
  //cmulti->Divide(3,2,0,0);

  /*
  hpedcheck->Draw();
  c1->Update();
  gPad->WaitPrimitive();
  c1->SaveAs("LowGain pedestal check.pdf");
  */
  

  char plot_title[50];
  int nevents = fChain->GetEntries();
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
  //  legend->SetNColumns(2);
  TImage *img = TImage::Create();


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
    for(int i=0; i<injevents_perdac; i++){
      ADC[i] = 0;
      n[i] = i+1;
    }

    //==================== Loop over the events ====================
    Crosstalk(52);

  for(int ev = 0; ev < nevents ; ++ev){  
    fChain -> GetEntry(ev); //== Get HITcollection from root file event 

    // Update the dac number and show plots for every event at the dac number 
    if(ev%injevents_perdac == 0){
      dac++; 
      sprintf(plot_title,"DAC:%d",dac_ctrl[dac]);
      TGraph *g = new TGraph(injevents_perdac,n,ADC);
      g->SetTitle(plot_title);
      //g->Draw("A*");
      //c1->Update();
      //gPad->WaitPrimitive();
      delete g;
      test=0;
    }
      
    for(int i = 0 ; i < NSCA ; ++i){
      TS[i] = HITS->rollposition[i];
    }
    
    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){
      
      H = HITS->Hits.at(hit);
      if(!H.CCorNC) continue; // remove unconnected channel

      if(H.formatCH == 31)
	{
      
	  for(int sca = 0; sca < NSCA; ++sca){
	    if( TS[sca] == ts ){	    
	      //H.SCA_hg[sca] -= avg_HG[H.chip][H.ch][sca]; // pedestal subtraction
	      //H.SCA_lg[sca] -= avg_LG[H.chip][H.ch][sca];

	
	      NoisyChannel_ADC_H[dac] = H.SCA_hg[sca];
	    }
	  }
	}
      

      for (int j=0; j<6; j++){
	
	  
	//	if(H.formatCH == cross_ch[j] ){
	if(H.formatCH == 42){
	  for(int sca = 0; sca < NSCA; ++sca){
	    if( TS[sca] == ts ){	    
	      //H.SCA_hg[sca] -= avg_HG[H.chip][H.ch][sca]; // pedestal subtraction
	      //H.SCA_lg[sca] -= avg_LG[H.chip][H.ch][sca];
	      h->Fill(H.SCA_hg[sca]);
	      cout << H.SCA_hg[sca] << endl;
	      //	  if(Cut(H.SCA_hg[sca], sigma_HG[H.chip][H.ch][sca])==1){
	      //ADC_H[dac]+=H.SCA_hg[sca];
	      // }
	      // if(Cut(H.SCA_lg[sca], sigma_LG[H.chip][H.ch][sca])==1){
	      //ADC_L[dac]+=H.SCA_lg[sca];
	      // cout << dac << " " << hit << " " << H.ch << " "<< sigma_LG[H.chip][H.ch][sca]<< " " << H.SCA_lg[sca] << endl;
	      //}
	      //h_TOTS->Fill(H.TOTS);
	      //h_TOTF->Fill(H.TOTF);
	      //h_TOAR->Fill(H.TOAR);
	      //h_TOAF->Fill(H.TOAF);
	    
	      //TOTS[dac] = H.TOTS;
	      Crosstalk_ADC_H[j][dac] = H.SCA_hg[sca];
	      Crosstalk_ADC_L[j][dac] = H.SCA_lg[sca];
	      Crosstalk_TOTS[j][dac] = H.TOTS;

	      
	      
	      //cout << H.TOTS <<" " << H.TOTF << endl <<
	      //  H.TOAR <<" " <<  H.TOAF << endl;
	    }
	  }
	}
      }
      test++;    
    }
  }
  
  
    for(int i=0; i<injevents; i++){
      ADC_H[i] /= injevents_perdac;
      ADC_L[i] /= injevents_perdac;
    }

    //==================== End of Loop ====================
    //
    //==================== Draw Plots ====================
    
    //g[ts] = new TGraph(injevents,dac_ctrl,ADC_H);
    TGraph* gh = new TGraph(injevents,dac_ctrl,ADC_H);
    TGraph* gl = new TGraph(injevents,dac_ctrl,ADC_L);
    TGraph* gTOT = new TGraph(injevents,dac_ctrl,TOTS);
    TGraph* ghlratio = new TGraph(injevents,ADC_L,ADC_H);
    TGraph* gltratio = new TGraph(injevents,TOTS,ADC_L);
    TGraph** gcross_h = new TGraph*[cross_num];
    TGraph** gcross_l = new TGraph*[cross_num];
    TGraph** gcross_TOTS = new TGraph*[cross_num];
    TGraph* gnoisy_h = new TGraph(injevents,dac_ctrl,NoisyChannel_ADC_H);
    
    
    h->Draw();
    c1->Update();
    gPad->WaitPrimitive();
    
    
    
    //  gh->GetYaxis()->SetRangeUser(0,200);
    //gl->GetXaxis()->SetRangeUser(xmin,xmax);

    // Draw Plots for each SCA 
    /*
      sprintf(plot_title,"High Gain");  
      g[ts]->SetTitle(plot_title);
      g[ts]->GetXaxis()->SetTitle("DAC");
      g[ts]->GetYaxis()->SetTitle("ADC");
      int style = 22+ts; int color = ts+1;
      if(color%10==0)color++;
      g[ts]->SetMarkerColor(color);
      g[ts]->SetMarkerStyle(style);
      //g[ts]->Draw("AP");
      //c1->Update();
      //gPad->WaitPrimitive();
      */
    //gPad->WaitPrimitive();

    
    /*h_TOTF->Draw();
      c1->Update();
      gPad->WaitPrimitive();
    
      h_TOAR->Draw();
      c1->Update();
      gPad->WaitPrimitive();
    
      h_TOAF->Draw();
      c1->Update();
      gPad->WaitPrimitive();
      h->Draw();
      c1->Update();
      gPad->WaitPrimitive();*/

    sprintf(plot_title,"Noisy Channel 31");
    gnoisy_h->SetTitle(plot_title);
    gnoisy_h->GetXaxis()->SetTitle("Event");
    gnoisy_h->GetYaxis()->SetTitle("ADC");
    gnoisy_h->SetMarkerStyle(7);
    gnoisy_h->Draw("AP");
    c1->Update();
    gPad->WaitPrimitive();

    sprintf(plot_title,"High Gain");
    gh->SetTitle(plot_title);
    gh->GetXaxis()->SetTitle("Event");
    gh->GetYaxis()->SetTitle("ADC");
    gh->SetMarkerStyle(7);
    gh->Draw("AP");
    c1->Update();
    //gPad->WaitPrimitive();
    
    sprintf(plot_title,"Low Gain");
    gl->SetTitle(plot_title);
    gl->GetXaxis()->SetTitle("Event");
    gl->GetYaxis()->SetTitle("ADC");
    gl->SetMarkerStyle(7);
    gl->Draw("AP");
    c1->Update();
    //gPad->WaitPrimitive();

    sprintf(plot_title,"TOTS");
    gTOT->SetTitle(plot_title);
    gTOT->GetXaxis()->SetTitle("Event");
    gTOT->GetYaxis()->SetTitle("ADC");
    gTOT->SetMarkerStyle(7);
    gTOT->Draw("AP");
    c1->Update();
    //gPad->WaitPrimitive();


    sprintf(plot_title,"High vs Low Gain");
    ghlratio->SetTitle(plot_title);
    ghlratio->GetXaxis()->SetTitle("Low");
    ghlratio->GetYaxis()->SetTitle("High");
    ghlratio->SetMarkerStyle(24);
    ghlratio->Draw("AP");
    c1->Update();
    //gPad->WaitPrimitive();

    sprintf(plot_title,"Low Gain vs TOT");
    gltratio->SetTitle(plot_title);
    gltratio->GetXaxis()->SetTitle("TOT");
    gltratio->GetYaxis()->SetTitle("Low");
    gltratio->SetMarkerStyle(24);
    gltratio->Draw("AP");
    c1->Update();
    //gPad->WaitPrimitive();

    for(int j=0; j<6; j++){
      if(cross_ch[j]!=257){
	for(int i=0; i<injevents; i++){
	  gcross_h[j] = new TGraph(injevents,dac_ctrl,Crosstalk_ADC_H[j]);
	  gcross_l[j] = new TGraph(injevents,dac_ctrl,Crosstalk_ADC_L[j]);
	  gcross_TOTS[j] = new TGraph(injevents,dac_ctrl,Crosstalk_TOTS[j]);
	}
	sprintf(plot_title,"HG_%d",cross_ch[j]);
	gcross_h[j]->SetTitle(plot_title);
	gcross_h[j]->SetMarkerStyle(20+j);
	gcross_h[j]->SetMarkerColor(2+j);
	gcross_h[j]->SetMarkerSize(0.5);
	gcross_h[j]->Draw("AP");
	c1->Update();
	gPad->WaitPrimitive();
	img->FromPad(c1);
	sprintf(img_title,"plots/crosstalk/InsertChannel=122_%s.png",plot_title);
	//img->WriteImage(img_title);

	sprintf(plot_title,"LG_%d",cross_ch[j]);
	gcross_l[j]->SetTitle(plot_title);
	gcross_l[j]->SetMarkerStyle(20+j);
	gcross_l[j]->SetMarkerColor(2+j);
	gcross_l[j]->SetMarkerSize(0.5);
	gcross_l[j]->Draw("AP");
	c1->Update();
	gPad->WaitPrimitive();
	img->FromPad(c1);
	sprintf(img_title,"plots/crosstalk/InsertChannel=123_%s.pdf",plot_title);
	//	img->WriteImage(img_title);
	c1->SaveAs(img_title);

	sprintf(plot_title,"TOTS_%d",cross_ch[j]);
	gcross_TOTS[j]->SetTitle(plot_title);
	gcross_TOTS[j]->SetMarkerStyle(20+j);
	gcross_TOTS[j]->SetMarkerColor(2+j);
	gcross_TOTS[j]->SetMarkerSize(0.5);
	gcross_TOTS[j]->Draw("AP");
	c1->Update();
	gPad->WaitPrimitive();
	img->FromPad(c1);
	sprintf(img_title,"plots/crosstalk/InsertChannel=122_%s.png",plot_title);
	//img->WriteImage(img_title);


	/*sprintf(plot_title,"LG_%d",cross_ch[j]);
	gl->SetTitle(plot_title);
	gl->Draw("AP");
	c1->Update();
	gPad->WaitPrimitive();
	
	sprintf(plot_title,"TOT_%d",cross_ch[j]);
	gTOT->SetTitle(plot_title);
	gTOT->Draw("AP");
	c1->Update();
	gPad->WaitPrimitive();*/

	mg->Add(gcross_h[j]);
	sprintf(leg,"HG_%d",cross_ch[j]);
	legend->AddEntry(gcross_h[j],leg,"p");	

      }
    }
    //    cmulti->Update();
    //cmulti->Draw();
    //gPad->WaitPrimitive();
    //img->FromPad(cmulti);
    //img->WriteImage("plots/chip3.png");
    mg->Draw("AP");
    gPad->Modified();
    legend->SetTextFont(110);
    legend->Draw();
    c1->Update();
    gPad->WaitPrimitive();
    c1->SetTitle("crosstalk");
    c1->SaveAs("crosstalk.pdf");
    img->FromPad(c1);
    //img->WriteImage("plots/crosstalk/InsertChannel=122_Crosstalk_low.png");



    /*
      mg->Add(g[ts]);
      sprintf(leg,"%dTS",ts);
      legend->AddEntry(g[ts],leg, "P");
    
      // h->Draw();
      //c1->Update();
      //gPad->WaitPrimitive();
      //}
      }
      //mg->SetTitle("High Gain");


      mg->Draw("AP");
      mg->GetXaxis()->SetTitle("DAC");
      mg->GetYaxis()->SetTitle("ADC");
      gPad->Modified();
      legend->Draw();
      c1->Update();
      //img->FromPad(c1);
      //img->WriteImage("MultiGraph0-6.png");
      gPad->WaitPrimitive();
    */
        
  //=================== End of filling hist =======================

    TwoDPoly();   
   
  // An example shows how to check the ADC vs TS for certain channel(injected)
  // Remove this comment to run

  
  TGraph *gr;
  for(int ev = 0; ev < nevents ; ++ev){
    if(ev % 35 != 0) continue;
    fChain -> GetEntry(ev);
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITS->rollposition[i];    
    int nhits = HITS->hit_num;
    for(int hit = 0; hit < nhits ; ++hit){
      H = HITS->Hits.at(hit);
      if(H.ch != 40) continue;

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
      //gPad->WaitPrimitive();
    }
  }
  
}

void makePlots::TwoDPoly(){
  
  //========== Here shows the example how to use TH2Poly to draw plot =========
  int nevents = fChain->GetEntries();
  char plot_title[50];
  TCanvas* c1 = new TCanvas();
  for(int ev = 0; ev < nevents ; ++ev){
    if(ev % 1 != 0) continue;

    fChain -> GetEntry(ev); // Get HITcollection from root file event
    TH2Poly *poly = new TH2Poly;
    InitTH2Poly(*poly); 
        
    for(int i = 0 ; i < NSCA ; ++i) TS[i] = HITS->rollposition[i];
    
    int nhits = HITS->hit_num;
    //    cout << nhits << endl;
    
    for(int hit = 0; hit < nhits ; ++hit){  
      H = HITS->Hits.at(hit);
      if(!H.CCorNC) continue; // remove unconnected channel
      int forCH = H.chip*32+H.ch/2;
      float X = CHmap[forCH].first;
      float Y = CHmap[forCH].second;	  

      if(H.ch==26){
	//if(H.ch==52 || forCH==65 || forCH==62 || forCH==63 || forCH ==61 || forCH==31 || forCH==95 || forCH==127){
	poly->Fill(X,Y,0);
      }
      else{
	for(int sca = 0; sca < NSCA; ++sca){
	  if(TS[sca] == 4 ){
	    H.SCA_hg[sca] -= avg_HG[H.chip][H.ch][sca]; // pedestal subtraction
	    H.SCA_hg[sca] += 50;
	    //H.SCA_lg[sca] -= avg_LG[H.chip][H.ch][sca];

	    poly->Fill(X,Y,H.SCA_hg[sca]);
	    //	    poly->Fill(X,Y,H.ch);
	  }
	}
      }
    }
    
    poly->Draw("colztext0");

    sprintf(plot_title,"HG_TS0_evt%d",ev);
    poly->SetTitle(plot_title);
    c1->Update();
    gPad->WaitPrimitive();
    delete poly;
  }
  //============== End of the example to use TH2Poly to draw plot ==============

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
  int sumHG     [NCHIP][NCH][NSCA];
  int sumLG     [NCHIP][NCH][NSCA];
  double sumsqHG   [NCHIP][NCH][NSCA];
  double sumsqLG   [NCHIP][NCH][NSCA];

  for(int i = 0; i < NCHIP; ++i){    
    for(int j = 0; j < NCH; ++j){
      for(int k = 0; k < NSCA; ++k){
	sumHG     [i][j][k] = 0;
	sumsqHG   [i][j][k] = 0;
	sumLG     [i][j][k] = 0;
	sumsqLG   [i][j][k] = 0;
	mem_of_SCA[i][j][k] = 0;}}}
  
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
      H = HITS->Hits.at(hit);
      for(int sca = 0; sca < NSCA; ++sca){
	if(TS[sca] >= skip_TS) continue;
	sumHG  [H.chip][H.ch][sca] += H.SCA_hg[sca];
	sumsqHG[H.chip][H.ch][sca] += H.SCA_hg[sca]*H.SCA_hg[sca];
	sumLG  [H.chip][H.ch][sca] += H.SCA_lg[sca];
	sumsqLG[H.chip][H.ch][sca] += H.SCA_lg[sca]*H.SCA_lg[sca];
	mem_of_SCA[H.chip][H.ch][sca]++;
      }
    }
  }

  for(int i = 0; i < NCHIP; ++i){    
    for(int j = 0; j < NCH; ++j){
      for(int k = 0; k < NSCA; ++k){
	avg_HG  [i][j][k] = (float)sumHG  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_HG[i][j][k] = (float)sumsqHG  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_HG[i][j][k] -= avg_HG  [i][j][k]*avg_HG  [i][j][k];
	sigma_HG[i][j][k] = sqrt(sigma_HG[i][j][k]);
	avg_LG  [i][j][k] = (float)sumLG  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_LG[i][j][k] = (float)sumsqLG  [i][j][k]/mem_of_SCA[i][j][k];
	sigma_LG[i][j][k] -= avg_LG  [i][j][k]*avg_LG  [i][j][k];
	sigma_LG[i][j][k] = sqrt(sigma_LG[i][j][k]); }}}
  
  if(output == true){
    string filename;
    filename = input_RUN;
    int start = filename.find_last_of("/");
    int end   = filename.find(".root");
    string outf = filename.substr(start+1,end-start-1);
    char outtitleH[100];
    char outtitleL[100];
    sprintf(outtitleH,"ped_result/%s_HG.txt",outf.c_str());
    ofstream fileHG(outtitleH);
    sprintf(outtitleL,"ped_result/%s_LG.txt",outf.c_str());
    ofstream fileLG(outtitleL);
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
	  fileHG << fixed << setprecision(2) << avg_HG[i][j][k] << " ";
	  fileLG << fixed << setprecision(2) << avg_LG[i][j][k] << " ";}
	fileHG << "\n";
	fileLG << "\n";
	fileHG << i << "\t" << j << "\t";
	fileLG << i << "\t" << j << "\t";

	for(int k = 0; k < NSCA; ++k){
	  fileHG << fixed << setprecision(2) << sigma_HG[i][j][k] << " ";
	  fileLG << fixed << setprecision(2) << sigma_LG[i][j][k] << " ";}
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
	  inHG >> avg_HG[chip][ch][sca];
	inHG >> chip >> ch;
	for(int sca = 0; sca < NSCA ; ++sca)
	  inHG >> sigma_HG[chip][ch][sca];
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
	  inLG >> avg_LG[chip][ch][sca];
	inLG >> chip >> ch;
	for(int sca = 0; sca < NSCA ; ++sca)
	  inLG >> sigma_LG[chip][ch][sca];
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
