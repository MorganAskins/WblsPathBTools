#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>
#include <TGraph.h>

#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/PathFit.hh>
#include <RAT/DS/Centroid.hh>
#include <RAT/DS/PMTInfo.hh>

using namespace std;

double dwall(TVector3&, double);
void marcfile(string, string, string, double, int);
TH2D MakeCDF(TH2D h, string name);

void marcfile(string iname, string oname, string hname, double radius, int subevchoice, int replace)
{
  // Input file
  unique_ptr<TFile> tfile(new TFile(iname.c_str()));
  TTree* t = (TTree*)tfile->Get("output");
  TTree* m = (TTree*)tfile->Get("meta");

  // Start with meta
  double eff = 0;
  m->SetBranchAddress("livetime", &eff);
  m->GetEvent(0);

  // Main Tree
  double x, y, z, chi2, mcke;
  int n100, n400, subev;
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);
  t->SetBranchAddress("mcke", &mcke);
  t->SetBranchAddress("n100", &n100);
  t->SetBranchAddress("chi2", &chi2);
  t->SetBranchAddress("n400", &n400);
  t->SetBranchAddress("subev", &subev);

  // replace?
  unique_ptr<TFile> gfile(new TFile("ibdgraph.root"));
  TGraph* hpool = (TGraph*)gfile->Get("hartlepool");
  TGraph* heysham = (TGraph*)gfile->Get("heysham");
  TGraph* bkg = (TGraph*)gfile->Get("background");

  double hpool_norm = hpool->Integral();
  double heysham_norm = heysham->Integral();
  double bkg_norm = bkg->Integral();

  // Histograms
  TH2D hist("pdf", "pdf", 90, 0, 9, 1000, 0, 1000);
  hist.Sumw2();
  int events_simulated = 0;
  int events_passed = 0;

  for(unsigned ev=0; ev<t->GetEntries(); ev++)
  {
    t->GetEvent(ev);
    if(subev == subevchoice || subev == -1)
    {
      events_simulated++;
      if(chi2>1.3) continue;
      events_passed++;
      TVector3 pos(x, y, z);
      double dw = dwall(pos, radius)/1000.0;
      // replace with heysham or bkg?
      double weight = 1.0;
      if( replace == 1 && subevchoice == 0) //heysham
        weight = heysham->Eval(mcke+1.8) / hpool->Eval(mcke+1.8) * hpool_norm / heysham_norm;
      if( replace == 2 && subevchoice == 0)
        weight = bkg->Eval(mcke+1.8) / hpool->Eval(mcke+1.8) * hpool_norm / bkg_norm;

      hist.Fill(dw, n400, weight);
    }
  }
  
  // New file
  unique_ptr<TFile> ofile( new TFile(oname.c_str(), "update") );

  TH2D cdf = MakeCDF(hist, hname);
  // Apply a scale factor.
  double simulated_volume = 2*pow((radius-150), 3)*TMath::Pi();
  double tank_volume = 2*pow(10000, 3)*TMath::Pi();
  double full_efficiency = eff * events_passed / events_simulated * simulated_volume / tank_volume;
  cdf.Scale(full_efficiency);

  ofile->Write(0, TObject::kOverwrite);
}

TH2D MakeCDF(TH2D h, string name)
{
  int xbins   = h.GetNbinsX();
  int ybins   = h.GetNbinsY();
  double xmin = h.GetXaxis()->GetBinUpEdge(0);
  double ymin = h.GetYaxis()->GetBinUpEdge(0);
  double xmax = xmin + h.GetXaxis()->GetBinWidth(0)*xbins;
  double ymax = ymin + h.GetYaxis()->GetBinWidth(0)*ybins;
  TH2D cdf(name.c_str(), name.c_str(), xbins, xmin, xmax, ybins, ymin, ymax);
  cdf.Sumw2();
  TH2D cdfy("cdf_y", "cdf_y", xbins, xmin, xmax, ybins, ymin, ymax);
  cdfy.Sumw2();
  double norm = h.GetSum();
  // reverse cdf integrate in one dimension
  for(int x=xbins-1; x>=0; x--)
  {
    double cumulative = 0;
    for(int y=ybins-1; y>=0; y--)
    {
      cumulative += h.GetBinContent(x, y) / norm;
      cdfy.SetBinContent(x, y, cumulative);
    }
  }
  // Then other dimension
  for(int y=ybins-1; y>=0; y--)
  {
    double cumulative = 0;
    for(int x=xbins-1; x>=0; x--)
    {
      cumulative += cdfy.GetBinContent(x, y);
      cdf.SetBinContent(x, y, cumulative);
    }
  }
  return cdf;
}

double dwall( TVector3 &v, double radius )
{
  // Arbitrary at the moment, could consider grabbing the pmt positions from
  // the run table.
  double top  =  radius;
  double bot  = -radius;
  double side =  radius;

  double rho = sqrt( v.X()*v.X() + v.Y()*v.Y() );
  double z = v.Z();
  double dtop  = top - z;
  double dbot  = z - bot;
  double dside = side - rho;

  double dcap = min(dtop, dbot);
  return min(dcap, dside);
}

int main(int argc, char** argv)
{
  // One file at a time, specify in / out
  // marcmap infile outfile radius subev
  if( argc == 6 )
    marcfile(argv[1], argv[2], argv[3], stod(argv[4]), stoi(argv[5]), 0);
  else
    marcfile(argv[1], argv[2], argv[3], stod(argv[4]), stoi(argv[5]), stoi(argv[6]));
  return 0;
}
