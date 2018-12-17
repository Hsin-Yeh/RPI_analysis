#ifndef PlotSetting_h
#define PlotSetting_h

#include "TROOT.h"
#include "TH2Poly.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <string>
using namespace std;

class PlotSetting{
 public:

  PlotSetting ();
  ~PlotSetting ();

  char plotfolder_path[200];

  //public function
  void GStd(TGraph& g, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot); 
  void G(TGraph& g, char* plot_title, string Xtitle, string Ytitle,int MarkerStyle,int MarkerColor,
		       int MarkerSize, int LineColor, int LineWidth, string Option,bool OptStat, bool Wait, bool SavePlot);
  void Multi(TMultiGraph& g, TLegend& legend, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot);
  void MultiAdd(TMultiGraph& multig, TGraph& g, TLegend& legend, char* plot_title, int MarkerStyle,int MarkerColor, int MarkerSize);
  void Poly(TH2Poly& poly, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot);
  void root_logon();
  int Color(int c);

 private:
  
};


#endif
