#include <vector>
#include <iostream>
#include <string>

#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>

#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Centroid.hh>
using namespace std;

void consolidate_fitter(string iname, string oname);

int main(int argc, char** argv)
{
  consolidate_fitter(argv[1], argv[2]);
  return 0;
}

void consolidate_fitter(string iname, string oname)
{
  TFile* tfile = new TFile(iname.c_str());
  TTree* T = (TTree*)tfile->Get("T");
  TTree* runT = (TTree*)tfile->Get("runT");

  RAT::DS::Run* run = new RAT::DS::Run();
  runT->SetBranchAddress("run", &run);

  RAT::DS::Root* ds = new RAT::DS::Root();
  T->SetBranchAddress("ds", &ds, 0);
  int entries = T->GetEntries();

  std::vector<double> *Q_Reco_X = new std::vector<double>;
  std::vector<double> *Q_Reco_Y = new std::vector<double>;
  std::vector<double> *Q_Reco_Z = new std::vector<double>;

  T->SetBranchAddress("Q_Reco_X", &Q_Reco_X);
  T->SetBranchAddress("Q_Reco_Y", &Q_Reco_Y);
  T->SetBranchAddress("Q_Reco_Z", &Q_Reco_Z);

  // Storage File / Tree
  TFile* otfile = new TFile(oname.c_str(), "recreate");
  TTree* newT = T->CloneTree(0);
  TTree* newR = runT->CloneTree(0);

  runT->GetEvent(0);
  newR->Fill();

  for(int i=0; i<entries; i++)
  {
    T->GetEvent(i);
    for(int j=0; j<ds->GetEVCount(); j++)
    {
      RAT::DS::EV* ev = ds->GetEV(j);
      RAT::DS::Centroid* centroid = ev->GetCentroid();
      centroid->SetPosFitName("ChargeFitter");
      TVector3 tv( Q_Reco_X->at(j), Q_Reco_Y->at(j), Q_Reco_Z->at(j) );
      centroid->SetPosition( tv );
    }

    newT->Fill();
  }

  otfile->Write(0, TObject::kOverwrite);
}
