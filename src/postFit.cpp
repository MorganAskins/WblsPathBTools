// Quick analysis of the data after running the fitting algorithm
// Mostly the same as postRat

// C++ Std
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;
// Boost
#include <boost/filesystem.hpp>
// Root
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
// Rat
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/PathFit.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCPhoton.hh>

int NhitsX(RAT::DS::EV* ev, double minT, double maxT );
vector<string> listDir(string directory);
string name_my_file(string input);
string last_dir(string input);
double dwall( TVector3 &v );

int main(int argc, char** argv)
{
  // Expectation is to receive a directory containing a particular
  // target with all generators nested. We loop over this directory.
  string dirname = "lowlevel_analysis/";
  string input_directory = argv[1];
  vector<string> target_directories = listDir(input_directory);

  // In our dirname, create one file, based on input_directory name
  string tfile_name = dirname + name_my_file(input_directory);

  RAT::DS::Root* ds = new RAT::DS::Root();

  TFile* outfile = new TFile( tfile_name.c_str(), "recreate" );

  for( auto td : target_directories )
  {
    TChain ch("T");
    ch.Add( (td+"/*").c_str() );
    ch.SetBranchAddress("ds", &ds);
    // Setup new tree
    string treename = last_dir(td);
    TTree* tree = new TTree(treename.c_str(), "pass 2");
    tree->SetDirectory( outfile );
    // Raw MC / Trigger
    TH1F* ped_hist = new TH1F( ("pedestal_"+treename).c_str(), ("pedestal_"+treename).c_str(), 1000, 0, 1000 );
    TH1F* mchits_hist= new TH1F( ("mchits_"+treename).c_str(), ("mchits_"+treename).c_str(), 1000, 0, 1000 );
    TH1F* n100_hist = new TH1F( ("n100_"+treename).c_str(), ("n100_"+treename).c_str(), 1000, 0, 1000 );
    TH1F* n400_hist = new TH1F( ("n400_"+treename).c_str(), ("n400_"+treename).c_str(), 1000, 0, 1000 );
    TH1F* hittimes_hist = new TH1F( ("hittimes_"+treename).c_str(), ("hittimes_"+treename).c_str(), 3000, -1000, 2000 );
    TH1F* trigs_hist = new TH1F( ("trigs_"+treename).c_str(), ("trigs_"+treename).c_str(), 10000, 0, 1000000 );
    TH1F* delta_hist = new TH1F( ("delta_"+treename).c_str(), ("delta_"+treename).c_str(), 10000, 0, 1000000 );
    TH2F* rho2z = new TH2F( ("rho2z_"+treename).c_str(), ("rho2z_"+treename).c_str(), 
        1000, 0, 10000*10000, 2000, -10000, 10000 );
    tree->Branch("ped_hist", &ped_hist);
    tree->Branch("mchits_hist", &mchits_hist);
    tree->Branch("n100_hist", &n100_hist);
    tree->Branch("n400_hist", &n400_hist);
    tree->Branch("hittimes_hist", &hittimes_hist);
    tree->Branch("trigs_hist", &trigs_hist);
    tree->Branch("delta_hist", &delta_hist);
    tree->Branch("rho2z", &rho2z);
    // Reconstruction
    TH1F* dx_hist = new TH1F( ("dx_"+treename).c_str(), ("dx_"+treename).c_str(), 4000, -2000, 2000 );
    TH1F* dy_hist = new TH1F( ("dy_"+treename).c_str(), ("dy_"+treename).c_str(), 4000, -2000, 2000 );
    TH1F* dz_hist = new TH1F( ("dz_"+treename).c_str(), ("dz_"+treename).c_str(), 4000, -2000, 2000 );
    TH2F* rec_rho2z = new TH2F( ("rec_rho2z_"+treename).c_str(), ("rec_rho2z_"+treename).c_str(), 
        1000, 0, 10000*10000, 2000, -10000, 10000 );
    TH2F* rec_xy = new TH2F( ("rec_xy"+treename).c_str(), ("rec_xy"+treename).c_str(), 
        2000, -10000, 10000, 2000, -10000, 10000 );
    TH1F* dwall_hist = new TH1F( ("dwall_"+treename).c_str(), ("dwall_"+treename).c_str(), 5000, -2000, 8000 );
    TH1F* dwall3_hist = new TH1F( ("dwall3_"+treename).c_str(), ("dwall3_"+treename).c_str(), 5000, 0, 2.0 );
    TH1F* good_hist = new TH1F(("good_"+treename).c_str(), ("good_"+treename).c_str(), 10000, 0, 1000);
    TH1F* dir_hist = new TH1F( ("dir_"+treename).c_str(), ("dir_"+treename).c_str(), 100, -1, 1 );
    TH2F* good2r = new TH2F( ("good2r_"+treename).c_str(), ("good2r_"+treename).c_str(), 
        1000, 0, 10000, 1000, 0, 500 );
    // Possibly add residuals

    tree->Branch("dx_hist", &dx_hist);
    tree->Branch("dy_hist", &dy_hist);
    tree->Branch("dz_hist", &dz_hist);
    tree->Branch("dwall_hist", &dwall_hist);
    tree->Branch("dwall3_hist", &dwall3_hist);
    tree->Branch("good_hist", &good_hist);
    tree->Branch("dir_hist", &dir_hist);
    tree->Branch("rec_rho2z", &rec_rho2z);
    tree->Branch("rec_xy", &rec_xy);
    tree->Branch("good2r", &good2r);
    for( int i=0; i<ch.GetEntries(); i++ )
    {
      ch.GetEvent(i);
      RAT::DS::MC* mc = ds->GetMC();
      int mchits = 0;
      for(int pmtid=0; pmtid < mc->GetMCPMTCount(); pmtid++)
      {
        RAT::DS::MCPMT* pmt = mc->GetMCPMT(pmtid);
        for(int photonid=0; photonid < pmt->GetMCPhotonCount(); photonid++)
        {
          RAT::DS::MCPhoton* photon = pmt->GetMCPhoton(photonid);
          double fectime = photon->GetFrontEndTime();
          hittimes_hist->Fill( fectime );
          if( !photon->IsDarkHit() ) mchits++;
        }
      }
      mchits_hist->Fill(mchits);
      TVector3 mcpos = mc->GetMCParticle(0)->GetPosition();
      TVector3 mcmom = mc->GetMCParticle(0)->GetMomentum();
      for(int subev=0; subev < ds->GetEVCount(); subev++)
      {
        RAT::DS::EV* ev = ds->GetEV(subev);
        RAT::DS::PathFit* pfit = ev->GetPathFit();
        TVector3 pos = pfit->GetPosition();
        TVector3 dir = pfit->GetDirection();

        dx_hist->Fill( mcpos.X() - pos.X() );
        dy_hist->Fill( mcpos.Y() - pos.Y() );
        dz_hist->Fill( mcpos.Z() - pos.Z() );
        double dR = sqrt( pow(mcpos.X()-pos.X(),2) + pow(mcpos.Y()-pos.Y(),2)* pow(mcpos.Z()-pos.Z(),2) );
        double dwall_calc = dwall( pos );
        dwall_hist->Fill( dwall_calc );
        dwall3_hist->Fill( pow((6700-dwall_calc), 3)/pow(6700,3) );

        double ddotd = (mcmom * dir) / mcmom.Mag() / dir.Mag();
        dir_hist->Fill(ddotd);

        double rho2 = pos.X()*pos.X() + pos.Y()*pos.Y();
        rec_rho2z->Fill( rho2, pos.Z() );
        rec_xy->Fill( pos.X(), pos.Y() );

        good_hist->Fill( pfit->GetGoodness() );
        good2r->Fill(dR, pfit->GetGoodness() );

        int n100 = NhitsX(ev, -20, 80);
        int n400 = NhitsX(ev, -20, 380);
        int ped = NhitsX(ev, -150, -50);
        ped_hist->Fill( ped );
        n100_hist->Fill( n100 );
        n400_hist->Fill( n400 );
        trigs_hist->Fill( ev->GetCalibratedTriggerTime() );
        delta_hist->Fill( ev->GetDeltaT() );
      }
      double mcrho2 = mcpos.X()*mcpos.X() + mcpos.Y()*mcpos.Y();
      double mcz = mcpos.Z();
      rho2z->Fill( mcrho2, mcz );
    }
    tree->Fill();
    cout << "Completed " << td << endl;
  }

  outfile->Write(0, TObject::kOverwrite);
  outfile->Close();
}

double dwall( TVector3 &v )
{
  // Arbitrary at the moment, could consider grabbing the pmt positions from
  // the run table.
  double top  =  6700;
  double bot  = -6700;
  double side =  6700;

  double rho = sqrt( v.X()*v.X() + v.Y()*v.Y() );
  double z = v.Z();
  double dtop  = top - z;
  double dbot  = z - bot;
  double dside = side - rho;

  double dcap = min(dtop, dbot);
  return min(dcap, dside);
}

int NhitsX(RAT::DS::EV* ev, double minT, double maxT )
{
  // Count hits relative to trigger between minT, and maxT
  int nhits = 0;
  for(int pmtc=0; pmtc < ev->GetPMTCount(); pmtc++)
  {
    double hit_time = ev->GetPMT(pmtc)->GetTime();
    if( hit_time > minT && hit_time < maxT )
      nhits++;
  }
  return nhits;
}

vector<string> listDir(string directory)
{
  vector<string> files;
  for(auto &p : boost::filesystem::directory_iterator( directory ) )
  {
    files.push_back( p.path().string() );
  }
  return files;
}

string name_my_file(string input)
{
  string nn = input;
  while( nn.find("/") != string::npos )
    nn.replace( nn.find("/"), 1, "_" );
  if( nn.back() == '_' )
    nn.replace( nn.length()-1, 1, ".root");
  else
    nn += ".root";
  return nn;
}

string last_dir(string input)
{
  string nn = input;
  while( nn.find("/") != string::npos )
    nn.replace(0, nn.find("/")+1, "");
  return nn;
}
