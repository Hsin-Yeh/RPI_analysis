// Arthor: Hsin-Yeh Wu
// Email : thankyouyou06@gmail.com


#ifndef makePlots_h
#define makePlots_h

#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH2Poly.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TExec.h"
#include "TStyle.h"
#include "PlotSetting.h"
#include <string>
#include <utility> //std::pair
#include <map>     //std::map

using namespace std;

const int NCHIP = 4;
const int NCH = 64;
const int NSCA = 13;
const int NformatCH = 128;
const int NCHANNEL = 256;
const int NRings = 5;

class makePlots{
 public:

  makePlots (TChain* inchain);
  ~makePlots();

  //public function
  void Init( string pedfile, string gainfile );
  void PlotProducer();
  void Pulse_display( int displayChannel = -1 , int acq_type = 0, int lowerR = -1, int upperR = -1, bool subPed_flag = true );

  //public parameter
  string input_fileName;
  int injCh;
  int ModuleNumber;
  double LG2HG_Conversion[NCHIP][NCH], TOT2LG_Conversion[NCHIP][NCH], HGTP[NCHIP][NCH], LGTP[NCHIP][NCH], TOTOffSet[NCHIP][NCH], ADC2MIP = 0.0227, LGTP_default = 900;
  PlotSetting P;

  
 private:
  
  void    yamlReader();
  void    GainFactorReader( string gainfile );
  double  mipConverter( double hg_SubPed, double lg_SubPed, double tot , int channel);
  int     ringPositionFinder( int inj_channel, int channel);
  double  CMCalculator( double **sig_SubPed, int *TS );
  bool    mipSigCheck( double *sig, int *TS );
  void    pulsePlotter( double *sig, int *TS, int ev, int ichip, int ich, int lowerR, int upperR );
  void    Gain_factor_producer();
  virtual Int_t    Cut(Long64_t entry, Long64_t sigma);
  void    Crosstalk(Int_t ch);
  void    P_and_N(int option,bool output);
  // P_and_N function:
  // option 0 can be used for pedestal run.
  // option 1 is an informal way that can deal with signal runs,
  // similar method is applied in test beam framework.
  double sigmaCal(int N, int sum, double sqsum);
  void   read_P_and_N(string ped_file);
  void   readmap();
  void   InitTH2Poly(TH2Poly& poly); //Give frame to TH2Poly
  

  TApplication   *app;
  TCanvas        *c;
  TTree          *Chain1;
  float          avg_HG[NCHIP][NCH][NSCA];
  float          sigma_HG[NCHIP][NCH][NSCA];
  float          avg_LG[NCHIP][NCH][NSCA];
  float          sigma_LG[NCHIP][NCH][NSCA];
  int            cross_ch_FirstRing[NCHIP][6];

  
  ///////////////////////////////
  // Declaration of leaf types //
  ///////////////////////////////

  Int_t         event;
  Int_t         chip;
  Int_t         roll;
  Int_t         dacinj;
  Int_t         timesamp[13];
  Int_t         hg[13][64];
  Int_t         lg[13][64];
  Int_t         tot_fast[64];
  Int_t         tot_slow[64];
  Int_t         toa_rise[64];
  Int_t         toa_fall[64];
  
  // map < key = chip*32+ch/2 , pair <x, y> > 
  map<int,pair < double,double > > CHmap;
};

#endif
