/////////////////////////////////////////////////////////
// Arthor: Chia-hung Chien  chchien521@gmail.com       
// Just use the same class name as we used to.         
/////////////////////////////////////////////////////////


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

class makePlots{
 public:

  makePlots (TChain* inchain);
  ~makePlots();

  //public function
  void Init();
  void PlotProducer();
  void Evt_display();
  void Inj_Pulse_display();
  void IdentifyInjCh();


  //public parameter
  bool Is_TB;
  string         input_RUN;
  bool doTruth;
  int pedopt;
  int Inj_ch;
  int ModuleNumber;
  double LG2HG_Conversion[NCHIP][NCH], TOT2LG_Conversion[NCHIP][NCH], HGTP[NCHIP][NCH], LGTP[NCHIP][NCH], TOTOffSet[NCHIP][NCH], ADC2MIP = 0.0227;
  double LGTP_default = 900;
  PlotSetting Plot;

  
 private:
  
  void yamlReader();
  void GainFactorReader();
  void Gain_factor_producer();
  virtual Int_t    Cut(Long64_t entry, Long64_t sigma);
  void Crosstalk(Int_t ch);
  void Crosstalk_2ndRing(Int_t ch);
  void P_and_N(int option,bool output);
  // P_and_N function:
  // option 0 can be used for pedestal run.
  // option 1 is an informal way that can deal with signal runs,
  // similar method is applied in test beam framework.
  double sigmaCal(int N, int sum, double sqsum);
  void read_P_and_N(string ped_file);
  void readmap();
  void InitTH2Poly(TH2Poly& poly); //Give frame to TH2Poly
  

  TFile* outfile;
  TApplication *app;
  TTree          *Chain1;
  int            TS[NSCA];
  float avg_HG_SCA  [NCHIP][NCH][NSCA];
  float sigma_HG_SCA[NCHIP][NCH][NSCA];
  float avg_LG_SCA  [NCHIP][NCH][NSCA];
  float sigma_LG_SCA[NCHIP][NCH][NSCA];
  int cross_ch_FirstRing[NCHIP][6];
  bool cross_type[NCHIP][6];
  int cross_ch_2ndRing[12];

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
