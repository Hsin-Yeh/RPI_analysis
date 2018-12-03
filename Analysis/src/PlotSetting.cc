#include "PlotSetting.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"

//Constructo
PlotSetting::PlotSetting()
{
}
//Destructor
PlotSetting::~PlotSetting()
{
}



void PlotSetting::TGraphPlotSetting(TGraph& g, char* plot_title, string Xtitle, string Ytitle,int MarkerStyle,int MarkerColor,
				  int MarkerSize, int LineColor, int LineWidth, string Option,bool OptStat, bool Wait, bool SavePlot) 
{
  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(OptStat);
  g.SetTitle(plot_title);
  g.GetXaxis()->SetTitle(Xtitle.c_str());
  g.GetYaxis()->SetTitle(Ytitle.c_str());
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
    sprintf(title,"%s/%s.pdf",plotfolder_path,plot_title);
    c->SaveAs(title);
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
  legend.Draw();
  c->Update();
  if(Wait){
    gPad->WaitPrimitive();
  }
  if(SavePlot){
    char title[200];
    sprintf(title,"%s/%s.pdf",plotfolder_path,plot_title);
    c->SaveAs(title);
    cout << plotfolder_path << endl;
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
  poly.Draw(Option.c_str());
  c->Update();
  if(Wait){
    gPad->WaitPrimitive();
  }
  if(SavePlot){
    char title[200];
    sprintf(title,"%s/%s.pdf",plotfolder_path,plot_title);
    c->SaveAs(title);
  }
  delete c;
}
