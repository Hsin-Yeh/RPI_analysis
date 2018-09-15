void long_term_ped(){
  
  TApplication* app;
  app = new TApplication("app",0,0);
  TCanvas* c1 = new TCanvas();
  char HG_name[100],LG_name[100];
  string filename;
  ifstream infile("long_term_ped.txt");
  int NCHIP = 4;
  int NCH = 256;
  int NSCA = 13;
  int count = 0;
  vector<double>  avg_HG, avg_LG, sigma_HG, sigma_LG, time;
  double  avg_HG_SCA[NCHIP][NCH][NSCA], avg_LG_SCA[NCHIP][NCH][NSCA], sigma_HG_SCA[NCHIP][NCH][NSCA], sigma_LG_SCA[NCHIP][NCH][NSCA];

  avg_HG.clear(); avg_LG.clear(); sigma_HG.clear(); sigma_LG.clear();

  //input pedestal text files 
  while(!infile.eof()){
    std::getline(infile, filename);
    if (filename.size()<2)
      {
	break;
      }
    filename.erase (filename.end()-7,filename.end());
    sprintf(HG_name,"%s_HG.txt",filename.c_str());
    sprintf(LG_name,"%s_LG.txt",filename.c_str());
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
      getline(inHG,line);
      avg_HG.push_back(atof(line.c_str()));
      getline(inHG,line);
      sigma_HG.push_back( atof(line.c_str()));

      getline(inLG,line);
      avg_LG.push_back(atof(line.c_str()));
      getline(inLG,line);
      sigma_LG.push_back(atof(line.c_str()));

      time.push_back(count);
      count++;
      /*
      inHG >> avg_HG;
      inHG >> sigma_HG;
      inLG >> avg_LG;
      inLG >> sigma_LG;
      */
      /*
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
	    inHG >> avg_HG_SCA[chip][ch][sca];
	  inHG >> chip >> ch;
	  for(int sca = 0; sca < NSCA ; ++sca)
	    inHG >> sigma_HG_SCA[chip][ch][sca];
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
	    inLG >> avg_LG_SCA[chip][ch][sca];
	  inLG >> chip >> ch;
	  for(int sca = 0; sca < NSCA ; ++sca)
	    inLG >> sigma_LG_SCA[chip][ch][sca];
	}      
      }
      */
      cout << "Reading pedestal file done!" << endl;
      inHG.close();
      inLG.close();
    }  
  }
  for (int i=0; i<avg_HG.size(); i++){
    cout << avg_HG.at(i) << " " << sigma_HG.at(i) << endl;
    cout << avg_LG.at(i) << " " << sigma_LG.at(i) << endl;
  }
  
  TGraphErrors* gavg_HG = new TGraphErrors(avg_HG.size(),&time[0],&avg_HG[0],0,&sigma_HG[0]);
  TGraph* gavg_LG = new TGraph(avg_HG.size(),&time[0],&avg_LG[0]);
  TGraph* gsigma_HG = new TGraph(avg_HG.size(),&time[0],&sigma_HG[0]);
  TGraph* gsigma_LG = new TGraph(avg_HG.size(),&time[0],&sigma_LG[0]);
  
  gavg_HG->SetMarkerStyle(3);
  gavg_HG->Draw("AP");
  c1->Update();
  gPad->WaitPrimitive();
  
  gavg_LG->SetMarkerStyle(3);
  gavg_LG->Draw("AP");
  c1->Update();
  gPad->WaitPrimitive();

  gsigma_HG->SetMarkerStyle(3);
  gsigma_HG->Draw("AP");
  c1->Update();
  gPad->WaitPrimitive();

  gsigma_LG->SetMarkerStyle(3);
  gsigma_LG->Draw("AP");
  c1->Update();
  gPad->WaitPrimitive();
  delete c1;
  
}
