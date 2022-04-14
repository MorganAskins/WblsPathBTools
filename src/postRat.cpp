// Quick analysis of the production data to check that everything
// is as expected. No reconstruction yet, just MC and trigger.

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
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCPhoton.hh>

int NhitsX(RAT::DS::EV* ev, double minT, double maxT);
vector<string> listDir(string directory);
string name_my_file(string input);
string last_dir(string input);

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
    TTree* tree = new TTree(treename.c_str(), "pass 1");
    tree->SetDirectory( outfile );
    TH1F* ped_hist = new TH1F( ("pedestal_"+treename).c_str(), ("pedestal_"+treename).c_str(), 1000, 0, 1000 );
    TH1F* n100_hist = new TH1F( ("n100_"+treename).c_str(), ("n100_"+treename).c_str(), 1000, 0, 1000 );
    TH1F* hittimes_hist = new TH1F( ("hittimes_"+treename).c_str(), ("hittimes_"+treename).c_str(), 3000, -1000, 2000 );
    TH1F* trigs_hist = new TH1F( ("trigs_"+treename).c_str(), ("trigs_"+treename).c_str(), 10000, 0, 1000000 );
    TH1F* delta_hist = new TH1F( ("delta_"+treename).c_str(), ("delta_"+treename).c_str(), 10000, 0, 1000000 );
    TH2F* rho2z = new TH2F( ("rho2z_"+treename).c_str(), ("rho2z_"+treename).c_str(), 
        1000, 0, 10000*10000, 2000, -10000, 10000 );
    tree->Branch("ped_hist", &ped_hist);
    tree->Branch("n100_hist", &n100_hist);
    tree->Branch("hittimes_hist", &hittimes_hist);
    tree->Branch("trigs_hist", &trigs_hist);
    tree->Branch("delta_hist", &delta_hist);
    tree->Branch("rho2z", &rho2z);
    for( int i=0; i<ch.GetEntries(); i++ )
    {
      ch.GetEvent(i);
      RAT::DS::MC* mc = ds->GetMC();
      for(int pmtid=0; pmtid < mc->GetMCPMTCount(); pmtid++)
      {
        RAT::DS::MCPMT* pmt = mc->GetMCPMT(pmtid);
        for(int photonid=0; photonid < pmt->GetMCPhotonCount(); photonid++)
        {
          RAT::DS::MCPhoton* photon = pmt->GetMCPhoton(photonid);
          double fectime = photon->GetFrontEndTime();
          hittimes_hist->Fill( fectime );
        }
      }
      for(int subev=0; subev < ds->GetEVCount(); subev++)
      {
        RAT::DS::EV* ev = ds->GetEV(subev);
        int n100 = NhitsX(ev, -20, 80);
        int ped = NhitsX(ev, -150, -50);
        ped_hist->Fill( ped );
        n100_hist->Fill( n100 );
        trigs_hist->Fill( ev->GetCalibratedTriggerTime() );
        delta_hist->Fill( ev->GetDeltaT() );
      }
      TVector3 mcpos = mc->GetMCParticle(0)->GetPosition();
      double rho2 = mcpos.X()*mcpos.X() + mcpos.Y()*mcpos.Y();
      double z = mcpos.Z();
      rho2z->Fill( rho2, z );
    }
    tree->Fill();
    cout << "Completed " << td << endl;
  }

  outfile->Write(0, TObject::kOverwrite);
  outfile->Close();
}

int NhitsX(RAT::DS::EV* ev, double minT, double maxT)
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
