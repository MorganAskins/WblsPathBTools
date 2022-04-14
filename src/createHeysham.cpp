#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom3.h>

#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PathFit.hh>

using namespace std;

void heysham(string iname, string oname, int option);

class graph
{
  public:
    TGraph* gr;
    vector<double> x;
    vector<double> y;
    graph(){}
    graph(TGraph* gr, vector<double> x)
    {
      this->gr = gr;
      this->x = x;
      for(auto &nx : x)
      {
        y.push_back(gr->Eval(nx));
      }
    }
    void scale(double value)
    {
      for(int i=0; i<y.size(); i++)
      {
        y[i] = y[i] * value;
      }
    }
    TGraph getTGraph()
    {
      TGraph grnew(x.size(), &x[0], &y[0]);
      return grnew;
    }
    graph divide(graph other)
    {
      vector<double> ny;
      for(int i=0; i<x.size(); i++)
      {
        ny.push_back(y[i] / other.y[i]);
      }
      graph newgraph;
      newgraph.x = x;
      newgraph.y = ny;
      return newgraph;
    }
};

int main(int argc, char** argv)
{
  heysham(argv[1], argv[2], stoi(argv[3]));
  return 0;
}

void heysham(string iname, string oname, int option)
{
  // Load the data
  TFile* tfile  = new TFile(iname.c_str());
  TTree* T      = (TTree*)tfile->Get("T");
  TTree* runT   = (TTree*)tfile->Get("runT");
  TTree* posdb  = (TTree*)tfile->Get("posdb");
  TTree* header = (TTree*)tfile->Get("header");

  RAT::DS::Run* run = new RAT::DS::Run();
  runT->SetBranchAddress("run", &run);

  int entries = T->GetEntries();
  RAT::DS::Root* ds = new RAT::DS::Root();
  T->SetBranchAddress("ds", &ds, 0);
  
  // TGraphs
  unique_ptr<TFile> gfile(new TFile("ibdgraph2.root"));
  TGraph* hpool = (TGraph*)gfile->Get("big_hartlepool");
  TGraph* convert;
  if( option == 0 )
    convert = (TGraph*)gfile->Get("heysham_signal");
  else
    convert = (TGraph*)gfile->Get("heysham_background");
  vector<double> xaxis;
  double min_x = 2.0;
  double max_x = 8.0;
  for(double i=2.0; i<8.0; i+=0.01)
    xaxis.push_back(i);

  graph hpoolgraph(hpool, xaxis);
  graph othergraph(convert, xaxis);
  // Since we are taking hartlepool and "oscillating", we cannot do better
  // in any given bin (stats wise) than that specific graph.
  double fraction = 100000;
  double bigdelta = -100000;
  double xchoice = 0;
  for(int i=0; i<xaxis.size(); i++)
  {
    double getfrac = hpoolgraph.y[i] / othergraph.y[i];
    double delta = othergraph.y[i] - hpoolgraph.y[i];
    if(getfrac < fraction)
    {
      fraction = getfrac;
      xchoice = xaxis[i];
      //printf("%f, %f\n", xchoice, fraction);
    }
  }
  printf("frac: %f @ %f\n", fraction, xchoice);
  othergraph.scale(fraction);

  // Storage File / Tree
  TFile* otfile = new TFile(oname.c_str(), "recreate");
  TTree* output = T->CloneTree(0);
  TTree* outrun = runT->CloneTree(0);
  TTree* outpos = posdb->CloneTree();
  TTree* outhead = header->CloneTree();

  TGraph testpool = hpoolgraph.getTGraph();
  testpool.SetName("pool");
  testpool.Write();
  TGraph testother = othergraph.getTGraph();
  testother.SetName("other");
  testother.Write();

  //graph oscillate = hpoolgraph.divide(othergraph);
  graph oscillate = othergraph.divide(hpoolgraph);
  TGraph osc = oscillate.getTGraph();
  osc.SetName("osc");
  osc.Write();

  runT->GetEvent(0);
  outrun->Fill();

  TRandom3* rndm = new TRandom3();

  // Loop through events
  for( int i=0; i < entries; i++ )
  {
    // Get New Event, if ANY sub event is within bounds do not cut.
    T->GetEvent(i);
    RAT::DS::MC* mc = ds->GetMC();
    // Find the positron if available
    int nparticles = mc->GetMCParticleCount();
    double energy = 0;
    for(int mcp=0; mcp<nparticles; mcp++)
    {
      RAT::DS::MCParticle* mpart = mc->GetMCParticle(mcp);
      int id = mpart->GetPDGCode();
      if( id == -11 )
      {
        energy = mpart->GetKE()+1.8;
      }
    }
    // Only apply to the first event.
    RAT::DS::EV* ev = ds->GetEV(0);
    double prob = osc.Eval(energy);
    if( energy < min_x || energy > max_x ) continue;
    if( rndm->Rndm() < prob ) output->Fill();
  }

  otfile->Write(0, TObject::kOverwrite);
}
