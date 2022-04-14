#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCPhoton.hh>

int main(int argc, char** argv)
{
  TFile* f = new TFile(argv[1]);
  TTree* t = (TTree*)f->Get("T");

  RAT::DS::Root* ds = new RAT::DS::Root();
  t->SetBranchAddress("ds", &ds);
  // Save
  TFile* out = new TFile("output_showtiming.root", "recreate");
  TH1F* hist = new TH1F("hist", "hist", 2000, 0, 2000);
  for(int i=0; i<t->GetEntries(); i++)
  {
    t->GetEvent(i);
    RAT::DS::MC* mc = ds->GetMC();
    for(int pmtid=0; pmtid < mc->GetMCPMTCount(); pmtid++)
    {
      RAT::DS::MCPMT* pmt = mc->GetMCPMT(pmtid);
      for(int photonid=0; photonid < pmt->GetMCPhotonCount(); photonid++)
      {
        RAT::DS::MCPhoton* photon = pmt->GetMCPhoton(photonid);
        double fectime = photon->GetFrontEndTime();
        hist->Fill(fectime);
      }
    }
  }
  out->Write(0, TObject::kOverwrite);
}
