#include "PlotSetting.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TImage.h"
#include "TMultiGraph.h"

//Constructo
PlotSetting::PlotSetting()
{
}
//Destructor
PlotSetting::~PlotSetting()
{
}


void PlotSetting::TGraphPlotStandard(TGraph& g, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot) 
{
  TCanvas* c = new TCanvas();
  g.SetTitle(plot_title);
  g.GetXaxis()->SetTitle(Xtitle.c_str());
  g.GetYaxis()->SetTitle(Ytitle.c_str());
  g.GetYaxis()->SetTitleOffset(1.3);
  g.Draw(Option.c_str());
  c->Update();
  if(Wait){
    gPad->WaitPrimitive();
  }
  if(SavePlot){
    char title[200];
    sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage(title);
  }
  delete c;  
}

void PlotSetting::TGraphPlotSetting(TGraph& g, char* plot_title, string Xtitle, string Ytitle,int MarkerStyle,int MarkerColor,
				    int MarkerSize, int LineColor, int LineWidth, string Option,bool OptStat, bool Wait, bool SavePlot) 
{
  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(OptStat);
  g.SetTitle(plot_title);
  g.GetXaxis()->SetTitle(Xtitle.c_str());
  g.GetYaxis()->SetTitle(Ytitle.c_str());
  g.GetYaxis()->SetTitleOffset(1.2);
  g.SetMarkerStyle(MarkerStyle);
  g.SetMarkerColor(MarkerColor);
  g.SetLineColor(LineColor);
  g.SetLineWidth(LineWidth);
  g.SetMarkerSize(MarkerSize);
  g.Draw(Option.c_str());
  c->Update();
  if(Wait){
    gPad->WaitPrimitive();
  }
  if(SavePlot){
    char title[200];
    sprintf(title,"%s/%s.png",plotfolder_path,plot_title);
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage(title);
  }
  delete c;  
}


void PlotSetting::TMultiGraphPlotSetting(TMultiGraph& g, TLegend& legend, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot) 
{
  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(OptStat);
  g.Draw(Option.c_str());
  g.SetTitle(plot_title);
  g.GetXaxis()->SetTitle(Xtitle.c_str());
  g.GetYaxis()->SetTitle(Ytitle.c_str());
  g.GetYaxis()->SetTitleOffset(1.4);
  legend.Draw();
  c->Update();
  if(Wait){
    gPad->WaitPrimitive();
  }
  if(SavePlot){
    char title[200];
    sprintf(title,"%s/%s.pdf",plotfolder_path,plot_title);
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage(title);
  }
  delete c;
}

void PlotSetting::TH2PolyPlotSetting(TH2Poly& poly, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot) 
{
  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(OptStat);
  poly.SetTitle(plot_title);
  poly.GetXaxis()->SetTitle(Xtitle.c_str());
  poly.GetYaxis()->SetTitle(Ytitle.c_str());
  poly.GetYaxis()->SetTitleOffset(1.4);
  poly.Draw(Option.c_str());
  c->Update();
  if(Wait){
    gPad->WaitPrimitive();
  }
  if(SavePlot){
    char title[200];
    sprintf(title,"%s/%s.pdf",plotfolder_path,plot_title);
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage(title);
  }
  delete c;
}

void PlotSetting::root_logon(){

  cout << endl << "Welcome to the ATLAS rootlogon.C" << endl;

  //
  // based on a style file from BaBar
  //

  //..BABAR style from RooLogon.C in workdir
  TStyle *atlasStyle= new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0;
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetStatColor(icol);
  //atlasStyle->SetFillColor(icol);

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);
  atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadTopMargin(0.1);
  atlasStyle->SetPadRightMargin(0.1);
  atlasStyle->SetPadBottomMargin(0.15);
  atlasStyle->SetPadLeftMargin(0.12);

  // use large fonts
  //Int_t font=72;
  Int_t font=42;
  Double_t tsize=0.045;
  atlasStyle->SetTextFont(font);



  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");

  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  atlasStyle->SetLabelOffset(tsize/4,"x");
  atlasStyle->SetLabelOffset(tsize/4,"y");
  


  //use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(0.5);
  atlasStyle->SetHistLineWidth(2.);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //get rid of X error bars and y error bar caps
  //atlasStyle->SetErrorX(0.001);

  //do not display any of the standard histogram decorations
  //atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  gROOT->SetStyle("Plain");

  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();



}

