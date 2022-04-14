#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TTimeStamp.h>
#include <TVector3.h>

#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/PathFit.hh>
#include <RAT/DS/Centroid.hh>
#include <RAT/DS/PMTInfo.hh>

using namespace std;

void ntuplefile(string iname, string oname);

int main(int argc, char** argv)
{
  // One file at a time, specify in / out
  ntuplefile(argv[1], argv[2]);
  return 0;
}

void ntuplefile(string iname, string oname)
{
  // Variables needed
  ULong64_t stonano = 1000000000;
  // Load the data
  TFile* tfile      = new TFile(iname.c_str());

  TTree* T          = (TTree*)tfile->Get("T");
  int entries       = T->GetEntries();
  RAT::DS::Root* ds = new RAT::DS::Root();
  string* dsname    = new string("signal");
  T->SetBranchAddress("ds", &ds, 0);
  if( T->GetListOfLeaves()->Contains("name") )
    T->SetBranchAddress("name", &dsname);

  // Storage File / Tree
  TFile* otfile = new TFile(oname.c_str(), "recreate");
  TTree* output = new TTree("output", "output");
  TTree* meta = new TTree("meta", "meta");

  // Branches to keep
  string name;
  double mcx, mcy, mcz;
  double mcu, mcv, mcw;
  double mcke;
  int evid;
  int subev;
  ULong64_t nanotime; // Should last 584 years
  // MCParticles
  int mcpcount;
  std::vector<Int_t> pdgcodes;
  std::vector<double> mcKEnergies;
  std::vector<double> mcPosx;
  std::vector<double> mcPosy;
  std::vector<double> mcPosz;
  std::vector<double> mcDirx;
  std::vector<double> mcDiry;
  std::vector<double> mcDirz;

  // Reconstructed variables / ev
  int pedestal;      // (-150, -50)
  int n100;          // (-20, 80)
  int n400;          // (-50, 350)
  int veto;          // (-50, 350) /type 2
  double Q;          // Integrated charge
  double maxQ;       // Maximum charge
  double vQ;         // Integrated charge veto
  double x, y, z;    // position
  double u, v, w;    // direction
  double chi2;       // goodness of fit
  double qx, qy, qz; // QFit
  int ncher;         // Cherenkov
  //int qvalid;     // QValid
  double qbal;
  double isotropyPath;
  double isotropyQfit;

  qx = 0;
  qy = 0;
  qz = 0;
  //qvalid = 1;

  output->Branch("name", &name);
  output->Branch("nanotime", &nanotime);
  output->Branch("mcx", &mcx);
  output->Branch("mcy", &mcy);
  output->Branch("mcz", &mcz);
  output->Branch("mcu", &mcu);
  output->Branch("mcv", &mcv);
  output->Branch("mcw", &mcw);
  output->Branch("mcke", &mcke);
  output->Branch("mcpcount", &mcpcount);
  output->Branch("evid", &evid);
  output->Branch("subev", &subev);
  output->Branch("pedestal", &pedestal);
  output->Branch("n100", &n100);
  output->Branch("n400", &n400);
  output->Branch("veto", &veto);
  output->Branch("Q", &Q);
  output->Branch("vQ", &vQ);
  output->Branch("maxQ", &maxQ);
  output->Branch("x", &x);
  output->Branch("y", &y);
  output->Branch("z", &z);
  output->Branch("u", &u);
  output->Branch("v", &v);
  output->Branch("w", &w);
  output->Branch("chi2", &chi2);
  output->Branch("ncher", &ncher);
  // QFit
  output->Branch("qx", &qx);
  output->Branch("qy", &qy);
  output->Branch("qz", &qz);
  // Classifiers
  output->Branch("qbal", &qbal);
  output->Branch("isotropyPath", &isotropyPath);
  output->Branch("isotropyQfit", &isotropyQfit);
  //output->Branch("qvalid", &qvalid);
  // Vector branches
  output->Branch("pdg", &pdgcodes);
  output->Branch("mcEnergy", &mcKEnergies);
  output->Branch("mcposx", &mcPosx);
  output->Branch("mcposy", &mcPosy);
  output->Branch("mcposz", &mcPosz);
  output->Branch("mcposx", &mcDirx);
  output->Branch("mcposy", &mcDiry);
  output->Branch("mcposz", &mcDirz);

  // Store for meta
  vector<double> vPed; // pedestals

  printf("Begin loop\n");
  // Loop through events
  int debug_line = 0;
  for( int i=0; i < entries; i++ )
  {
    // Clear old event
    pdgcodes.clear();
    mcKEnergies.clear();
    mcPosx.clear();
    mcPosy.clear();
    mcPosz.clear();
    mcDirx.clear();
    mcDiry.clear();
    mcDirz.clear();
    // Get New Event
    T->GetEvent(i);
    RAT::DS::MC* mc = ds->GetMC();
    TTimeStamp mcTTS = mc->GetUTC();
    ULong64_t mctime = static_cast<ULong64_t>(mcTTS.GetSec())*stonano +
                       static_cast<ULong64_t>(mcTTS.GetNanoSec());
    name = *dsname;
    mcpcount = mc->GetMCParticleCount();
    // Get MC Particle Information
    for( int p=0; p<mcpcount; p++ )
    {
      RAT::DS::MCParticle* particle = mc->GetMCParticle(p);
      pdgcodes.push_back( particle->GetPDGCode() );
      mcKEnergies.push_back( particle->GetKE() );
      TVector3 mcpos = particle->GetPosition();
      TVector3 mcdir = particle->GetMomentum();
      mcPosx.push_back( mcpos.X() );
      mcPosy.push_back( mcpos.Y() );
      mcPosz.push_back( mcpos.Z() );
      mcDirx.push_back( mcdir.X()/mcdir.Mag() );
      mcDiry.push_back( mcdir.Y()/mcdir.Mag() );
      mcDirz.push_back( mcdir.Z()/mcdir.Mag() );
    }
    mcx = mcPosx[0];
    mcy = mcPosy[0];
    mcz = mcPosz[0];
    //mcu = mcDirx[0];
    //mcv = mcDiry[0];
    //mcw = mcDirz[0];
    mcke = accumulate(mcKEnergies.begin(), mcKEnergies.end(), 0.0);
    // Store aggregate particle info (first position, sum of ke)
    // Get Sub Events and write to ttree
    for( subev=0; subev < ds->GetEVCount(); subev++ )
    {
      RAT::DS::EV* ev = ds->GetEV(subev);
      evid = ev->GetID();
      nanotime = static_cast<ULong64_t>(ev->GetCalibratedTriggerTime()) + mctime;
      RAT::DS::PathFit* fit = ev->GetPathFit();
      TVector3 pos = fit->GetPosition();
      x = pos.X();
      y = pos.Y();
      z = pos.Z();
      TVector3 dir = fit->GetDirection();
      u = dir.X();
      v = dir.Y();
      w = dir.Z();
      chi2 = fit->GetGoodness();
      pedestal = 0;
      n100     = ev->GetTotalCharge();
      n400     = ev->GetTotalCharge();
      veto     = 0;
      Q        = ev->GetDeltaT(); 
      vQ       = ev->GetTotalCharge(); 
      maxQ     = ev->GetTotalCharge(); 
      ncher    = ev->GetDeltaT();
      vPed.push_back(pedestal);
      // QFit
      RAT::DS::Centroid* qfit = ev->GetCentroid();
      TVector3 qpos = qfit->GetPosition();
      qx = qpos.X();
      qy = qpos.Y();
      qz = qpos.Z();

      // Fill
      output->Fill();
    }
  }
  //
  double avg_pedestal = std::accumulate( vPed.begin(), vPed.end(), 0.0 ) / vPed.size();
  vector<double> diff(vPed.size());
  transform(vPed.begin(), vPed.end(), diff.begin(),
      [avg_pedestal](double x){return x - avg_pedestal;});
  double sq_sum = inner_product( diff.begin(), diff.end(), diff.begin(), 0.0 );
  double std_pedestal = sqrt(sq_sum / vPed.size());
  printf("Avg %f, std %f \n", avg_pedestal, std_pedestal);
  meta->Branch("AvgPedestal", &avg_pedestal);
  meta->Branch("StdPedestal", &std_pedestal);
  // If header, store livetime, else make one up
  double livetime = -1.0;;
  if( tfile->GetListOfKeys()->Contains("header") )
  {
    TTree* header = (TTree*)tfile->Get("header");
    // Now which type of "header"?
    if( header->GetListOfLeaves()->Contains("livetime") )
    {
      // This is a merged background
      header->SetBranchAddress("livetime", &livetime);
    }
    else if( header->GetListOfLeaves()->Contains("efficiency") )
    {
      header->SetBranchAddress("efficiency", &livetime);
    }
    header->GetEvent(0);
  }
  meta->Branch("livetime", &livetime);
  //
  meta->Fill();

  otfile->Write(0, TObject::kOverwrite);
}
