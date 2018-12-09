#ifndef PlotSetting_h
#define PlotSetting_h

#include "TROOT.h"
#include "TH2Poly.h"
#include "TLegend.h"
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
  void TGraphPlotStandard(TGraph& g, char* plot_title, string Xtitle, string Ytitle, string Option, bool Wait, bool SavePlot); 
  void TGraphPlotSetting(TGraph& g, char* plot_title, string Xtitle, string Ytitle,int MarkerStyle,int MarkerColor,
		       int MarkerSize, int LineColor, int LineWidth, string Option,bool OptStat, bool Wait, bool SavePlot);
  void TMultiGraphPlotSetting(TMultiGraph& g, TLegend& legend, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot);
  void TH2PolyPlotSetting(TH2Poly& poly, char* plot_title, string Xtitle, string Ytitle, string Option,bool OptStat, bool Wait, bool SavePlot);
  void root_logon();

 private:
  
};


#endif
