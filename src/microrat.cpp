#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>
#include <MicroDS.hh>

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

using namespace std;

// This is a compromise between the ntuple and ratds data structure
// which will preserve some of the properties of the ratds (so not flat
// or easy to parse) but reduced in size.

void microfile(string iname, string oname);
int NhitsX(RAT::DS::EV* ev, double minT, double maxT);

int main(int argc, char** argv)
{
  // One file at a time, specify in / out
  microfile(argv[1], argv[2]);
  return 0;
}

void microfile(string iname, string oname)
{
  // Load the data
  TFile* tfile = new TFile(iname.c_str());
  TTree* T = (TTree*)tfile->Get("T");
  TTree* head = (TTree*)tfile->Get("header");

  int entries = T->GetEntries();
  RAT::DS::Root* ds = new RAT::DS::Root();
  string* dsname = new string();
  T->SetBranchAddress("ds", &ds, 0);
  // T->SetBranchAddress("name", &dsname);

  // Storage File / Tree
  TFile* otfile = new TFile(oname.c_str(), "recreate");
  TTree* output = new TTree("micro", "micro");
  TTree* meta   = new TTree("meta", "meta");
  TTree* header = head->CloneTree(0);
  head->GetEvent(0);
  header->Fill();

  MicroDS* mds = new MicroDS();
  mds->NewBranches( output );

  // Store for meta
  vector<double> vPed; // pedestals

  // Loop through events
  for( int i=0; i < entries; i++ )
  {
    mds->cleardata();

    // Get New Event
    T->GetEvent(i);
    RAT::DS::MC* mc = ds->GetMC();
    mds->mcT = mc->GetUTC();
    //name = *dsname;
    mds->mcpcount = mc->GetMCParticleCount();
    // Get MC Particle Information
    for( int p=0; p<mds->mcpcount; p++ )
    {
      RAT::DS::MCParticle* particle = mc->GetMCParticle(p);
      mds->pdgcodes->push_back( particle->GetPDGCode() );
      mds->mcKEnergies->push_back( particle->GetKE() );
      TVector3 mcpos = particle->GetPosition();
      TVector3 mcdir = particle->GetMomentum();
      mds->mcPosx->push_back( mcpos.X() );
      mds->mcPosy->push_back( mcpos.Y() );
      mds->mcPosz->push_back( mcpos.Z() );
      mds->mcDirx->push_back( mcdir.X()/mcdir.Mag() );
      mds->mcDiry->push_back( mcdir.Y()/mcdir.Mag() );
      mds->mcDirz->push_back( mcdir.Z()/mcdir.Mag() );
    }
    // Store aggregate particle info (first position, sum of ke)
    // Get Sub Events and write to ttree
    for( int sub=0; sub < ds->GetEVCount(); sub++ )
    {
      RAT::DS::EV* ev = ds->GetEV(sub);
      mds->subev->push_back(sub);
      mds->triggertime->push_back( ev->GetCalibratedTriggerTime() );
      RAT::DS::PathFit* fit = ev->GetPathFit();
      TVector3 pos = fit->GetPosition();
      mds->x->push_back(pos.X());
      mds->y->push_back(pos.Y());
      mds->z->push_back(pos.Z());
      TVector3 dir = fit->GetDirection();
      mds->u->push_back(dir.X());
      mds->v->push_back(dir.Y());
      mds->w->push_back(dir.Z());
      int pedcount = NhitsX(ev, -150, -50);
      mds->pedestal->push_back(pedcount);
      mds->n100->push_back( NhitsX(ev, -20, 80) );
      mds->n400->push_back( NhitsX(ev, -50, 350) );
      vPed.push_back(pedcount);
      // Fill
    }
    output->Fill();
  }
  //
  double avg_pedestal = std::accumulate( vPed.begin(), vPed.end(), 0.0 ) / vPed.size();
  double std_pedestal = sqrt( std::inner_product( vPed.begin(), vPed.end(), vPed.end(), 0.0 ) / vPed.size() -
      avg_pedestal*avg_pedestal );
  meta->Branch("AvgPedestal", &avg_pedestal);
  meta->Branch("StdPedestal", &std_pedestal);
  meta->Fill();

  otfile->Write(0, TObject::kOverwrite);
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
