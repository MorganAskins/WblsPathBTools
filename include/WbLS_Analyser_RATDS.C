///------------ Performs Q_Fit on RAT-PAC output -------------////////
///------------ Author: Leon Pickard             -------------////////
///------------ Email: ljpickard@ucdavis.edu     -------------////////
///------------ Date: 25/06/2020                 -------------////////
///------------ To run the macro, you should provide the root filename and state whether you wish to produce a PDF, or a fit using the PDF.root file (e.g. .x WbLS_Analyser_RATDS.C("output.root",1) to perform the fit on output.root). Note, wildcard entries can be used when forming the PDF.  

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMTInfo.hh>
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---- Function to perform the minimisation. This method voxelises the detector (it isn't particularly fast...). To perform this minimisation the fitting argument is 1. ----//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
TVector3 Fitting_Likelihood_Volumes(RAT::DS::Root* ds, RAT::DS::PMTInfo* pmtinfo, RAT::DS::EV *ev, TH2F *PDF){

  //---- Reconstructed position vector
  TVector3 Position = {-1e9,-1e9,-1e9};
  
  if (ds->GetEVCount()){
    
    //---- Define some variables. Note, Jump, Iterations, Volume_Reduction_Factor and Jump_Scaling can be optimised for speed and precision. 
    float_t Likelihood_Best = -1e10; float_t Jump = 5000.0;  float_t Max_Range = 12000.0; float_t Iterations = 5; float_t Volume_Reduction_Factor = 1.5; float_t Jump_Scaling = 1.75;
    TVector3 Test_Vertex; TVector3 Best_Vertex; TVector3 Start_Vector = {-Jump,-Jump,-Jump}; TVector3 End_Vector = {Jump, Jump, Jump};
    
    //---- Scanning over the detector volume
    for (int Attempt=0; Attempt<Iterations; Attempt += 1){
      for (float_t x_pos=Start_Vector[0]; x_pos<End_Vector[0]+1; x_pos += Jump) {
	for (float_t y_pos=Start_Vector[1]; y_pos<End_Vector[1]+1; y_pos += Jump) {
	  if (sqrt(x_pos*x_pos+y_pos*y_pos) > Max_Range){continue;}
	  for (float_t z_pos=Start_Vector[2]; z_pos<End_Vector[2]+1; z_pos += Jump) {
	    Test_Vertex = {x_pos,y_pos,z_pos};	  
	    float_t Likelihood = 0;
	    
	    for(long iPMT = 0; iPMT < ev->GetPMTCount(); iPMT++ ){
	      int PMT_ID             = ev->GetPMT(iPMT)->GetID();
              int pmt_type           = pmtinfo->GetType(PMT_ID);
              if( pmt_type != 1 ) continue;
	      TVector3 PMT_Position  = pmtinfo->GetPosition(PMT_ID);
	      TVector3 PMT_Direction = pmtinfo->GetDirection(PMT_ID);	      
	      TVector3 R_Test_Vector = Test_Vertex - PMT_Position;
	      float_t Angle          = cos(R_Test_Vector.Angle(PMT_Direction));	      
	      Likelihood += ev->GetPMT(iPMT)->GetCharge()*log(PDF->GetBinContent(PDF->FindBin(R_Test_Vector.Mag(),Angle)));    
	    }    
	  
	    for (long ipmt = 0; ipmt < pmtinfo->GetPMTCount(); ipmt++){
              int pmt_type           = pmtinfo->GetType(ipmt);
              if( pmt_type != 1 ) continue;
	      TVector3 PMT_Position  = pmtinfo->GetPosition(ipmt);
	      TVector3 PMT_Direction = pmtinfo->GetDirection(ipmt);	      
	      TVector3 R_Test_Vector = Test_Vertex - PMT_Position;
	      float_t Angle          = cos(R_Test_Vector.Angle(PMT_Direction));
	      Likelihood -= PDF->GetBinContent(PDF->FindBin(R_Test_Vector.Mag(),Angle));
	    }

	   //---- If we find a test vertex with a larger likelihood, that is the new reconstructed position 
	  if (Likelihood > Likelihood_Best){ Likelihood_Best = Likelihood; Position = Test_Vertex;}
	  }
	}
      }

      //---- After each scan of the detector, refocus the scan near the best fit position and scan with a finer step size 
      for (int iAxis = 0; iAxis < 3; iAxis += 1){
	Start_Vector[iAxis] = Position[iAxis] - Jump/Volume_Reduction_Factor;
	End_Vector[iAxis]   = Position[iAxis] + Jump/Volume_Reduction_Factor;
      }
      Jump = Jump/Jump_Scaling;
    }
  }

  //---- After scanning for the number of iterations defined, the position with the largest likelihood gives the reconstructed position
  return Position;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------//
//---- Ascent function to perform minimisation. This method uses a random walk approach. To perform this minimisation the fitting argument is 2. ----//
//---------------------------------------------------------------------------------------------------------------------------------------------------//
TVector3 Fitting_Likelihood_Ascent(RAT::DS::Root* ds, RAT::DS::PMTInfo* pmtinfo, RAT::DS::EV *ev ,TH2F *PDF){

  //---- Reconstructed position vector
  TVector3 Position = {-1e9,-1e9,-1e9};
  TRandom3* gRandom = new TRandom3();
  
  if (ds->GetEVCount()){
    
    //---- Define some variables. 
    //float_t Likelihood_Best = -1e10; float Jump = 1000.0; int iWalk_Max = 70; TVector3 Test_Vertex = {0,0,0}; TVector3 Best_Vertex; float Start = 3000.0; TVector3 End_Vector = {0,0,0}; vector<TVector3> Vector_List = {{0,0,0},{Start,0,0},{-Start,0,0},{0,Start,0},{0,-Start,0},{0,0,Start},{0,0,-Start},{Start,Start,0},{Start,-Start,0},{-Start,Start,0},{-Start,-Start,0},{Start,Start,Start},{Start,-Start,Start},{-Start,Start,Start},{-Start,-Start,Start},{Start,Start,-Start},{Start,-Start,-Start},{-Start,Start,-Start},{-Start,-Start,-Start}};
    float_t Likelihood_Best = -1e10; float Jump = 2000.0; int iWalk_Max = 200; TVector3 Test_Vertex = {0,0,0}; TVector3 Best_Vertex; float Start = 5000.0; TVector3 End_Vector = {0,0,0}; vector<TVector3> Vector_List = {{0,0,0},{Start,0,0},{-Start,0,0},{0,Start,0},{0,-Start,0},{0,0,Start},{0,0,-Start},{Start,Start,0},{Start,-Start,0},{-Start,Start,0},{-Start,-Start,0},{Start,Start,Start},{Start,-Start,Start},{-Start,Start,Start},{-Start,-Start,Start},{Start,Start,-Start},{Start,-Start,-Start},{-Start,Start,-Start},{-Start,-Start,-Start}};
    
    int iWalkNoMore = 19;
    //int iWalkNoMore = 30;
    //---- Go for a walk...
    for(int iWalk = 0; iWalk < iWalk_Max; iWalk++){
      if (iWalk < iWalkNoMore){Test_Vertex = Vector_List[iWalk];}
      else if (iWalk == iWalkNoMore){gRandom->Sphere(End_Vector[0],End_Vector[1],End_Vector[2],Jump); Test_Vertex = Position + End_Vector;}

      float_t Likelihood = 0;
      
      for(long iPMT = 0; iPMT < ev->GetPMTCount(); iPMT++ ){
	int PMT_ID             = ev->GetPMT(iPMT)->GetID();
        int pmt_type           = pmtinfo->GetType(PMT_ID);
        if( pmt_type != 1 ) continue;
	TVector3 PMT_Position  = pmtinfo->GetPosition(PMT_ID);
	TVector3 PMT_Direction = pmtinfo->GetDirection(PMT_ID);	      
	TVector3 R_Test_Vector = Test_Vertex - PMT_Position;
	float_t Angle          = cos(R_Test_Vector.Angle(PMT_Direction));	      
	Likelihood += ev->GetPMT(iPMT)->GetCharge()*log(PDF->GetBinContent(PDF->FindBin(R_Test_Vector.Mag(),Angle)));    
      }    
      
      for (long ipmt = 0; ipmt < pmtinfo->GetPMTCount(); ipmt++){
        int pmt_type           = pmtinfo->GetType(ipmt);
        if( pmt_type != 1 ) continue;
	TVector3 PMT_Position  = pmtinfo->GetPosition(ipmt);
	TVector3 PMT_Direction = pmtinfo->GetDirection(ipmt);	      
	TVector3 R_Test_Vector = Test_Vertex - PMT_Position;
	float_t Angle          = cos(R_Test_Vector.Angle(PMT_Direction));
	Likelihood -= PDF->GetBinContent(PDF->FindBin(R_Test_Vector.Mag(),Angle));
      }

      //---- If we find a test vertex with a larger likelihood, that is the new reconstructed position 
      if (Likelihood > Likelihood_Best){Likelihood_Best = Likelihood; iWalk--; Jump=Jump/1.05; Position = Test_Vertex;
	if (End_Vector[0] != 0 && End_Vector[1] !=0 && End_Vector[2] !=0){Test_Vertex = Position + End_Vector;}
	else {gRandom->Sphere(End_Vector[0],End_Vector[1],End_Vector[2],Jump); Test_Vertex = Position + End_Vector;}
      }
      else{gRandom->Sphere(End_Vector[0],End_Vector[1],End_Vector[2],Jump); Test_Vertex = Position + End_Vector;}
    }
  }

  //---- After scanning for the number of iterations defined, the position with the largest likelihood gives the reconstructed position
  return Position;
}


//-----------------------------------------------------------------------------------------------------------------//
//---- Function to produce the Q fit PDF                                                                       ----//
//-----------------------------------------------------------------------------------------------------------------//
void PDF_Producer(const char* filename_ratpac){

  RAT::DSReader *dsReader;
  dsReader = new RAT::DSReader(filename_ratpac);
  
  TChain* runtri;
  runtri = new TChain("runT");
  
  if (TString(filename_ratpac).MaybeWildcard()) {
    //---- Assume there is a runT in all files
    runtri->Add(filename_ratpac);
    RAT::DS::RunStore::SetReadTree(runtri);
  }
  else { 	//---- In single file case, we can check
    TFile *ftemp = TFile::Open(filename_ratpac);
    if (ftemp->Get("runT")) {
      runtri->Add(filename_ratpac);
      RAT::DS::RunStore::SetReadTree(runtri);} //---- else, no runT, so don't register runtri with RunStore
    delete ftemp;}
  
  RAT::DS::Root *ds;
  RAT::DS::Run* run;

  ULong64_t NbEntries = dsReader->GetTotal();
  int Minimum_NHits   = 5;
  int X_NBins         = 300;
  int Y_NBins         = 200;
  int X_Min           = 0;
  int X_Max           = 30000;
  int Y_Min           = -1;
  int Y_Max           = 1;
  
  //---- Output PDF file
  TFile f_output("PDF_Output.root","RECREATE");
  
  //---- Create some histograms
  TH2F* h_R_Cos_Theta      = new TH2F("h_R_Cos_Theta","R_Cos_Theta",X_NBins,X_Min,X_Max,Y_NBins,Y_Min,Y_Max);
  TH2F* h_R_Cos_Theta_Hits = new TH2F("h_R_Cos_Theta_Hits","R_Cos_Theta_Hits",X_NBins,X_Min,X_Max,Y_NBins,Y_Min,Y_Max);

  //---- Analysis loop over all the events to produce PDF
  for (ULong64_t entry=0; entry<NbEntries; ++entry) {

    if (entry%100 == 0){cout << "Entry:    " << entry << " of " << NbEntries << endl;}
    
    ds  = dsReader->GetEvent(entry);
    run = RAT::DS::RunStore::Get()->GetRun(ds);
    RAT::DS::PMTInfo* pmtinfo = run->GetPMTInfo();

    TVector3 Interaction_Vertex = ds->GetMC()->GetMCParticle(0)->GetPosition();
    
    for(long iPMT = 0; iPMT < ds->GetMC()->GetMCPMTCount(); iPMT++ ){
      for(long iPhot = 0; iPhot < ds->GetMC()->GetMCPMT(iPMT)->GetMCPhotonCount(); iPhot++){
	int PMT_ID             = ds->GetMC()->GetMCPMT(iPMT)->GetID();
        int pmt_type           = pmtinfo->GetType(PMT_ID);
        if( pmt_type != 1 ) continue;
	TVector3 PMT_Position  = pmtinfo->GetPosition(PMT_ID);
	TVector3 PMT_Direction = pmtinfo->GetDirection(PMT_ID);	  
	TVector3 R_Vector      = Interaction_Vertex - PMT_Position;
	float_t Angle          = cos(R_Vector.Angle(PMT_Direction));
	
	if (ds->GetMC()->GetMCPMTCount() >= Minimum_NHits){
	  h_R_Cos_Theta->Fill(R_Vector.Mag(),Angle);
	}
      }
    }
    
    for (long ipmt = 0; ipmt < pmtinfo->GetPMTCount(); ipmt++){
      int pmt_type           = pmtinfo->GetType(ipmt);
      if( pmt_type != 1 ) continue;
      TVector3 PMT_Position  = pmtinfo->GetPosition(ipmt);
      TVector3 PMT_Direction = pmtinfo->GetDirection(ipmt);
      TVector3 R_Vector      = Interaction_Vertex - PMT_Position;
      float_t Angle          = cos(R_Vector.Angle(PMT_Direction));
      if (ds->GetMC()->GetMCPMTCount() > Minimum_NHits){
	h_R_Cos_Theta_Hits->Fill(R_Vector.Mag(),Angle);
      }  
    }
  }
  
  //---- Write PDF
  f_output.cd();
  h_R_Cos_Theta->Divide(h_R_Cos_Theta_Hits);
  h_R_Cos_Theta->Write();
  f_output.Close();
  
  return;
}

//-----------------------------------------------------------------------------------------------------------------//
//---- Function to add the Q Fit vertices to the exisitng RATDS file -----------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
void RATDS_Output(const char* filename_ratpac,int Fitting){

  RAT::DSReader *dsReader;
  dsReader = new RAT::DSReader(filename_ratpac);
  
  RAT::DS::EV *ev; 
  RAT::DS::Root *ds;
  RAT::DS::Run* run;
  
  TChain* runtri;
  runtri = new TChain("runT");
  
  runtri->Add(filename_ratpac);
  RAT::DS::RunStore::SetReadTree(runtri);

  ULong64_t NbEntries = dsReader->GetTotal();
  int Minimum_NHits   = 5;
  std::vector<double> *Q_Reco_X = new vector<double>;
  std::vector<double> *Q_Reco_Y = new vector<double>;
  std::vector<double> *Q_Reco_Z = new vector<double>;
  std::vector<int> *Q_Fit_Valid = new vector<int>;

  //---- Load the Q fit PDF
  TFile *f      = new TFile("PDF.root");
  TH2F * PDF    = (TH2F*)f->Get("h_R_Cos_Theta");

  //---- Update the RATDS output file by adding the Q_Fit result
  TFile *g           = new TFile(filename_ratpac,"update");
  TTree *T           = (TTree*)g->Get("T");
  
  //---- Check if the QFitter has run already
  TBranch *bQRX, *bQRY, *bQRZ, *Q_F_Valid;; 
  if( T->GetListOfLeaves()->Contains("Q_Fit_Valid") )
  {
    T->SetBranchAddress("Q_Reco_X", &Q_Reco_X, &bQRX);
    T->SetBranchAddress("Q_Reco_Y", &Q_Reco_Y, &bQRY);
    T->SetBranchAddress("Q_Reco_Z", &Q_Reco_Z, &bQRZ);
    T->SetBranchAddress("Q_Fit_Valid", &Q_Fit_Valid, &Q_F_Valid);
  }
  else
  {
    bQRX      = T->Branch("Q_Reco_X",&Q_Reco_X);
    bQRY      = T->Branch("Q_Reco_Y",&Q_Reco_Y);
    bQRZ      = T->Branch("Q_Reco_Z",&Q_Reco_Z);
    Q_F_Valid = T->Branch("Q_Fit_Valid",&Q_Fit_Valid);
  }
  
  //---- Analysis loop over all the events
  for (ULong64_t entry=0; entry<NbEntries; ++entry) {
    
    if (entry%10 == 0){cout << "Entry:    " << entry << " of " << NbEntries << endl;}
    
    ds  = dsReader->GetEvent(entry);
    run = RAT::DS::RunStore::Get()->GetRun(ds);
    RAT::DS::PMTInfo* pmtinfo = run->GetPMTInfo();

    Q_Reco_X->clear(); Q_Reco_Y->clear(); Q_Reco_Z->clear(); Q_Fit_Valid->clear();
    
    for (int iE=0; iE < ds->GetEVCount(); iE++){
      
      ev  = ds->GetEV(iE);
      
      if (ev->GetPMTCount() < Minimum_NHits) {Q_Reco_X->push_back(-1e9); Q_Reco_Y->push_back(-1e9); Q_Reco_Z->push_back(-1e9); Q_Fit_Valid->push_back(0); continue;}

      TVector3 Best_Fit;
      //---- Perform vertex fitting and fill the branches
      if (Fitting == 1){Best_Fit = Fitting_Likelihood_Volumes(ds, pmtinfo, ev, PDF);}
      else if (Fitting == 2){Best_Fit = Fitting_Likelihood_Ascent(ds, pmtinfo, ev, PDF);}
	
      Q_Reco_X->push_back(Best_Fit[0]);
      Q_Reco_Y->push_back(Best_Fit[1]);
      Q_Reco_Z->push_back(Best_Fit[2]);
      Q_Fit_Valid->push_back(1);      
    }

    bQRX->Fill();
    bQRY->Fill();
    bQRZ->Fill();
    Q_F_Valid->Fill();

    Q_Reco_X->clear(); Q_Reco_Y->clear(); Q_Reco_Z->clear(); Q_Fit_Valid->clear();
    
  }

  //---- Write reco positions to RATDS file and clean up
  T->Write("",TObject::kOverwrite); f->Close(); g->Close(); delete dsReader;//, run, tri, runtri;
  
  delete Q_Reco_X;
  delete Q_Reco_Y;
  delete Q_Reco_Z;
  delete Q_Fit_Valid;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------//
//---- Main function. First argument is(are) the filename(s), the second is whether you wish to produce a PDF from them or to perform the fitting. ----//
//-----------------------------------------------------------------------------------------------------------------------------------------------------//
void WbLS_Analyser_RATDS(const char* filename_ratpac, int Fitting) {
  
  //---- Load RAT libraries (for dsReader)
  // gSystem->Load("$(RATPAC_PATH)/lib/libRATEvent.so");
  
  //---- Initialisation
  // Clashes with the time.h included in rat
  // std::clock_t start;
  // double duration;

  //---- Starts the timer
  // start = clock();

  //---- Fit PDF or produce PDF?
  if (!Fitting){PDF_Producer(filename_ratpac);}
  else if (Fitting){RATDS_Output(filename_ratpac, Fitting);}
    
  //---- Ends the timer
  // duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
  // cout << "Execution time: " << duration << " seconds\n";

}
