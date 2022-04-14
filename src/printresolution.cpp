#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TTimeStamp.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TF1.h>

using namespace std;
void doprint(string iname);

int main(int argc, char** argv)
{
  // One file at a time, specify in / out
  doprint(argv[1]);
  return 0;
}

void doprint(string iname)
{
  unique_ptr<TFile> tfile(new TFile(iname.c_str()));
  TTree* t = (TTree*)tfile->Get("output");
  double x, y, z, mcx, mcy, mcz, mcke;
  int subev;
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);
  t->SetBranchAddress("mcx", &mcx);
  t->SetBranchAddress("mcy", &mcy);
  t->SetBranchAddress("mcz", &mcz);
  t->SetBranchAddress("mcke", &mcke);
  t->SetBranchAddress("subev", &subev);

  TH1D dx("dx", "dx", 200, -5000, 5000);
  TH1D dy("dy", "dy", 200, -5000, 5000);
  TH1D dz("dz", "dz", 200, -5000, 5000);

  unique_ptr<TF1> f(new TF1("Fittt", "[0]*TMath::Exp(-(x-[1])^2/2/[2]^2)", -5000, 5000));

  for(int i=0; i<t->GetEntries(); i++)
  {
    t->GetEvent(i);
    if( abs(mcx) > 4000 || abs(mcy) > 4000 || abs(mcz) > 4000 ) continue;
    if( subev != 0 ) continue;
    if( mcke < 2.0 ) continue;
    dx.Fill(x - mcx);
    dy.Fill(y - mcy);
    dz.Fill(z - mcz);
  }
  f->SetParameters(t->GetEntries()/100, 0, 100);
  dx.Fit(f.get(), "Q");
  double xv = abs(f->GetParameter(2));
  f->SetParameters(t->GetEntries()/100, 0, 100);
  dy.Fit(f.get(), "Q");
  double yv = abs(f->GetParameter(2));
  f->SetParameters(t->GetEntries()/100, 0, 100);
  dz.Fit(f.get(), "Q");
  double zv = abs(f->GetParameter(2));
  printf("%s %0.2f %0.2f %0.2f\n", iname.c_str(), xv, yv, zv);
}
