/////////////////////////////////////////////////////////
// Arthor: Chia-hung Chien  chchien521@gmail.com       
// Just use the same class name as we used to.         
/////////////////////////////////////////////////////////


#ifndef makePlots_h
#define makePlots_h

#include "hit_hits_class.h"
#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH2Poly.h"
#include "TApplication.h"
#include <string>
#include <utility> //std::pair
#include <map>     //std::map

using namespace std;
// Please Avoid variable NSCA = 13 and NCHIP = 4 defined in hit_hits_class.h
const int NCH = 64;
const int NformatCH = 128;

class makePlots{
 public:

  makePlots (TChain* inchain);
  ~makePlots();

  void Init();
  void PlotProducer();
  void Ntuplizer();
  void Evt_display();
  void Inj_Pulse_display();
  void IdentifyInjCh();
  string         input_RUN;
  bool doTruth;
  int pedopt;
  int Inj_ch;
  
 private:

  virtual Int_t    Cut(Long64_t entry, Long64_t sigma);
  void Crosstalk(Int_t ch);
  void Crosstalk_2ndRing(Int_t ch);
  void P_and_N(int option,bool output);

  // P_and_N function:
  // option 0 can be used for pedestal run.
  // option 1 is an informal way that can deal with signal runs,
  // similar method is applied in test beam framework.
  void Pedestal_ana(int option);
  void Injection_ana(int Inj_ch);
  double sigmaCal(int N, int sum, double sqsum);
  
  void read_P_and_N(string ped_file);
  void readmap();
  void InitTH2Poly(TH2Poly& poly); //Give frame to TH2Poly

  TFile* outfile;
  TApplication *app;
  TTree          *Chain1;
  hitcollection  *HITCOLLECTION;
  hit            H;
  int            TS[NSCA];
  float avg_HG_SCA  [NCHIP][NCH][NSCA];
  float sigma_HG_SCA[NCHIP][NCH][NSCA];
  float avg_LG_SCA  [NCHIP][NCH][NSCA];
  float sigma_LG_SCA[NCHIP][NCH][NSCA];
  int cross_ch[6];
  int cross_ch_2ndRing[12];
  std::vector<int> InjCh = {2};

  // map < key = chip*32+ch/2 , pair <x, y> > 
  map<int,pair < double,double > > CHmap;
};

#endif
