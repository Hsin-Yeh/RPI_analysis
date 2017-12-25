#include "hit_hits_class.h"

void hit::Fill_hit(bool _CCorNC,int _chtype,int _chip,int _ch,double _posx,double _posy,int _TOTS,int _TOTF,int _TOAR,int _TOAF,int _SCA_hg[13],int _SCA_lg[13]){

  CCorNC = _CCorNC;
  chtype = _chtype;
  chip   = _chip;
  ch     = _ch;
  formatCH = (_CCorNC) ? (chip*32+ch/2) : -1;
  TOTS = _TOTS;
  TOTF = _TOTF;
  TOAR = _TOAR;  
  TOAF = _TOAF;
  posx = _posx;
  posy = _posy;
  
  for(int i = 0 ; i < NSCA; ++i){
  SCA_hg[i] = _SCA_hg[i];
  SCA_lg[i] = _SCA_lg[i];}

}
