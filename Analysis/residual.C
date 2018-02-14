void residual(){
  char fname[100];
  int tarCH[8] = {0,2,10,20,30,40,50,60};
  int chs[5] = {1,2,4,6,8};
  float dac_a[1000];
  float HG_a[5][8][4][1000];
  float LG_a[5][8][4][1000];
  float TOT_a[5][8][4][1000];

  float dac;
  float HG[4],LG[4],TOT[4];
  
  for(int Id = 0; Id < 8 ; ++Id){

    for(int num = 0 ; num < 5 ; ++num){
      sprintf(fname,"Ntuple_CH%d_inj%d.root",tarCH[Id],chs[num]);
      if(Id == 1 && num == 3) continue;
      if(Id == 4 && num == 3) continue;
      
      TFile *f  = new TFile(fname);
      TTree *tt = (TTree *) f->Get("tree");

      tt->SetBranchAddress("dac",&dac);
      tt->SetBranchAddress("HG",HG);
      tt->SetBranchAddress("LG",LG);
      tt->SetBranchAddress("TOT",TOT);
      int nev = tt->GetEntries();
      
      for(int ev = 0; ev < nev ; ++ev){
	tt->GetEntry(ev);
	dac_a[ev] = dac;
	for(int ski = 0; ski < 4 ; ++ski){
	  HG_a[num][Id][ski][ev] = HG[ski];
	  LG_a[num][Id][ski][ev] = LG[ski];
	  TOT_a[num][Id][ski][ev]= TOT[ski];
	}  
      }
      f->Close();
    }
  }
  TCanvas *c1 = new TCanvas();
  TGraph *gr;
  
  char desc[100];
  vector<float> res;
  for(int Id = 0; Id < 8; ++Id){
    for(int Gain = 0; Gain < 3; ++Gain){
      for(int chip = 0 ; chip < 4 ; ++chip){

	TMultiGraph *mgr_h = new TMultiGraph();
	TLegend *leg = new TLegend(0.5,0.13,0.87,0.3);
	leg->SetBorderSize(0);
      
	for(int num = 1 ; num < 5 ; ++num){
	  if(Id == 1 && num == 3) continue;
	  if(Id == 4 && num == 3) continue;
	  res.clear();
	  for(int ev = 0; ev < 1000 ; ++ev){
	    float a = (HG_a[num][Id][chip][ev] - HG_a[0][Id][chip][ev] )/HG_a[0][Id][chip][ev];
	    float b = (LG_a[num][Id][chip][ev] - LG_a[0][Id][chip][ev] )/LG_a[0][Id][chip][ev];
	    float c = (TOT_a[num][Id][chip][ev]- TOT_a[0][Id][chip][ev] )/TOT_a[0][Id][chip][ev];
	    if(Gain == 0) res.push_back(a - 0.25*(num-2));
	    if(Gain == 1) res.push_back(b - 0.25*(num-2));
	    if(Gain == 2) res.push_back(c - 0.25*(num-2));
	  
	  }
	  gr = new TGraph(1000, dac_a, &res[0] );
	  gr->SetMarkerStyle(19+num);
	  gr->SetMarkerSize(0.8);
	  gr->SetMarkerColor(num);
	  mgr_h->Add(gr,"P");
	  sprintf(desc,"(inj%d-inj1)/inj1 + (%.2f)",chs[num],0.25*(num-2));
	  leg->AddEntry(gr,desc,"P");

	  /*
	    gr = new TGraph(1000, dac_a, &res_lg[0] );
	    mgr_l->Add(gr);
	    gr = new TGraph(1000, dac_a, &res_tot[0] );
	    mgr_t->Add(gr);
	  */

	  //base line
	  double X[2] = {0,4000};
	  double Y[2] = {(-0.25)*(num-2),(-0.25)*(num-2)};
	  TGraph *grl = new TGraph(2,X,Y);
	  grl->SetLineColor(num);
	  grl->SetLineWidth(3.5);
	  grl->SetLineStyle(7);
	  mgr_h->Add(grl,"L");
	}
	mgr_h->Draw("APL");
	mgr_h->GetXaxis()->SetTitle("dac");
	char pT[40];
	if(Gain == 0)
	  sprintf(pT,"HG_CH%d_chip%d_residual",tarCH[Id],chip);
	if(Gain == 1)
	  sprintf(pT,"LG_CH%d_chip%d_residual",tarCH[Id],chip);
	if(Gain == 2)
	  sprintf(pT,"TOT_CH%d_chip%d_residual",tarCH[Id],chip);
	mgr_h->GetYaxis()->SetTitle(pT);
	mgr_h->SetMaximum(0.5);
	mgr_h->SetMinimum(-1);
	leg->Draw("same");

	c1->Update();
	sprintf(pT,"%s.png",pT);
	c1->SaveAs(pT);
      } 
    }
  }


  //Non residual
    for(int Id = 0; Id < 8; ++Id){
    for(int Gain = 0; Gain < 3; ++Gain){
      for(int chip = 0 ; chip < 4 ; ++chip){

	TMultiGraph *mgr_h = new TMultiGraph();
	if(Gain <= 1)
	  TLegend *leg = new TLegend(0.5,0.13,0.87,0.4);
	else
	  TLegend *leg = new TLegend(0.5,0.63,0.87,0.87);
	leg->SetBorderSize(0);
      
	for(int num = 0 ; num < 5 ; ++num){
	  if(Id == 1 && num == 3) continue;
	  if(Id == 4 && num == 3) continue;

	  if(Gain == 0)
	    gr = new TGraph(1000, dac_a, HG_a[num][Id][chip] );
	  if(Gain == 1)
	    gr = new TGraph(1000, dac_a, LG_a[num][Id][chip] );
	  if(Gain == 2)
	    gr = new TGraph(1000, dac_a, TOT_a[num][Id][chip] );

	  gr->SetMarkerStyle(20+num);
	  gr->SetMarkerSize(0.8);
	  gr->SetMarkerColor(num+1);
	  if(num == 4)
	    gr->SetMarkerColor(num+2);

	  mgr_h->Add(gr,"P");
	  if(Gain == 0)
	    sprintf(desc,"HG_inj%d_CH%d_chip%d",chs[num],tarCH[Id],chip);
	  if(Gain == 1)
	    sprintf(desc,"LG_inj%d_CH%d_chip%d",chs[num],tarCH[Id],chip);
	  if(Gain == 2)
	    sprintf(desc,"TOT_inj%d_CH%d_chip%d",chs[num],tarCH[Id],chip);

	  leg->AddEntry(gr,desc,"P");
	}
	mgr_h->Draw("AP");
	mgr_h->GetXaxis()->SetTitle("dac");
	
	
	sprintf(pT,"ADC");
	mgr_h->GetYaxis()->SetTitle(pT);
	if(Gain == 2)
	  mgr_h->SetMaximum(1000);
	//mgr_h->SetMinimum(-1);
	leg->Draw("same");

	c1->Update();
	if(Gain == 0)
	  sprintf(pT,"HG_CH%d_chip%d.png",tarCH[Id],chip);
	if(Gain == 1)
	  sprintf(pT,"LG_CH%d_chip%d.png",tarCH[Id],chip);
	if(Gain == 2)
	  sprintf(pT,"TOT_CH%d_chip%d.png",tarCH[Id],chip);

	c1->SaveAs(pT);
	//getchar();
      } 
    }
  }

}
