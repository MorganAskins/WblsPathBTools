#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TVector3.h>

#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PathFit.hh>

using namespace std;

void prunefile(string iname, string oname, double x, double y, double z);
int NhitsX(RAT::DS::EV* ev, double minT, double maxT);
double dwall( TVector3 &v, double wallradius );
double dwall( TVector3 &v, double wallradius, double wallz );
double dwall( TVector3 &v, double wallx, double wally, double wallz );

// Assume we provide x y z
int main(int argc, char** argv)
{
  double x = stod(argv[3]);
  double y = stod(argv[4]);
  double z = stod(argv[5]);
  prunefile(argv[1], argv[2], x, y ,z);
  return 0;
}

void prunefile(string iname, string oname, double wallx, double wally, double wallz)
{
  // Load the data
  TFile* tfile = new TFile(iname.c_str());
  TTree* T = (TTree*)tfile->Get("T");
  TTree* runT = (TTree*)tfile->Get("runT");

  RAT::DS::Run* run = new RAT::DS::Run();
  runT->SetBranchAddress("run", &run);

  int entries = T->GetEntries();
  RAT::DS::Root* ds = new RAT::DS::Root();
  T->SetBranchAddress("ds", &ds, 0);

  // Storage File / Tree
  TFile* otfile = new TFile(oname.c_str(), "recreate");
  TTree* output = T->CloneTree(0);
  TTree* outrun = runT->CloneTree(0);
  TH1I* pedestal = new TH1I("pedestal", "pedestal", 1000, 0, 1000);

  runT->GetEvent(0);
  outrun->Fill();

  // Create a simple position database
  vector<double> xdb, ydb, zdb;
  TTree* dbtree = new TTree("posdb", "posdb");
  dbtree->Branch("xdb", &xdb);
  dbtree->Branch("ydb", &ydb);
  dbtree->Branch("zdb", &zdb);

  int keep_entries = 0;
  // Loop through events
  for( int i=0; i < entries; i++ )
  {
    // Get New Event, if ANY sub event is within bounds do not cut.
    T->GetEvent(i);
    double maxwall = -1000000;
    int max_n100 = 0;
    double max_goodness = 0;

    // since the posdb only matters for single, just keep the first
    double xfirst = -99999;
    double yfirst = -99999;
    double zfirst = -99999;
    for(int evc=0; evc<ds->GetEVCount(); evc++)
    {
      // Swapping to charge fitter (Centroid)
      RAT::DS::EV* ev = ds->GetEV(evc);
      RAT::DS::PathFit* pf = ev->GetPathFit();
      TVector3 pos = pf->GetPosition();
      double goodness = pf->GetGoodness();
      if(goodness > max_goodness) max_goodness = goodness;

      //RAT::DS::Centroid* centroid = ev->GetCentroid();
      //TVector3 pos = centroid->GetPosition();
      double dw = dwall(pos, wallx, wally, wallz);
      if( dw > maxwall ) maxwall = dw;
      //TVector3 dir = pf->GetDirection();
      int n100 = NhitsX(ev, -20, 80);
      if(n100 > max_n100) max_n100 = n100;
      int ped = NhitsX(ev, -150, -50);
      pedestal->Fill(ped);

      // pdb
      if( evc == 0 )
      {
        xfirst = pos.X();
        yfirst = pos.Y();
        zfirst = pos.Z();
      }
    }
    // Fiducial volume
    if( maxwall < 0 ) continue;
    // Nhits
    if( max_n100 < 5 ) continue;
    // Good fit
    if( max_goodness > 0.99 ) continue;
    xdb.push_back( xfirst );
    ydb.push_back( yfirst );
    zdb.push_back( zfirst );
    output->Fill();
    keep_entries++;
  }
  dbtree->Fill();

  TTree* header = new TTree("header", "header");
  double efficiency = static_cast<double>(keep_entries) / entries;
  header->Branch("efficiency", &efficiency);
  // Prune can be run multiple times, even on hadd'd files, to account for
  // this we check if the header tree already exists, if so we get the average
  // efficiency from there, otherwise we calculate a new one.
  if( tfile->GetListOfKeys()->Contains("header") )
  {
    TTree* oldheader = (TTree*)tfile->Get("header");
    efficiency = 0;
    double old_eff;
    oldheader->SetBranchAddress("efficiency", &old_eff);
    for(int i=0; i<oldheader->GetEntries(); i++)
    {
      oldheader->GetEvent(i);
      efficiency += old_eff / oldheader->GetEntries();
    }
  }
  header->Fill();

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

double dwall( TVector3 &v, double wallradius )
{
  // Arbitrary at the moment, could consider grabbing the pmt positions from
  // the run table.
  double top  =  wallradius;
  double bot  = -wallradius;
  double side =  wallradius;

  double rho = sqrt( v.X()*v.X() + v.Y()*v.Y() );
  double z = v.Z();
  double dtop  = top - z;
  double dbot  = z - bot;
  double dside = side - rho;

  double dcap = min(dtop, dbot);
  return min(dcap, dside);
}

double dwall( TVector3 &v, double wallradius, double wallz )
{
  // Arbitrary at the moment, could consider grabbing the pmt positions from
  // the run table.
  double top  =  wallz;
  double bot  = -wallz;
  double side =  wallradius;

  double rho = sqrt( v.X()*v.X() + v.Y()*v.Y() );
  double z = v.Z();
  double dtop  = top - z;
  double dbot  = z - bot;
  double dside = side - rho;

  double dcap = min(dtop, dbot);
  return min(dcap, dside);
}

double dwall( TVector3 &v, double wallx, double wally, double wallz )
{
  // Nearest wall
  double dx = wallx - abs(v.X());
  double dy = wally - abs(v.Y());
  double dz = wallz - abs(v.Z());
  return min(dx,min(dy,dz));
}
