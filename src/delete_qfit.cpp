#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <vector>
#include <string>
using namespace std;

int main(int argc, char** argv)
{
  vector<string> args(argv+1, argv+argc);
  // Fix files
  for(auto v : args)
  {
    TFile* f = new TFile(v.c_str(), "update");
    TTree* T = (TTree*)f->Get("T");
    TBranch *qx = T->GetBranch("Q_Reco_X");
    TBranch *qy = T->GetBranch("Q_Reco_Y");
    TBranch *qz = T->GetBranch("Q_Reco_Z");
    TBranch *qv = T->GetBranch("Q_Fit_Valid");
    T->GetListOfBranches()->Remove(qx);
    T->GetListOfBranches()->Remove(qy);
    T->GetListOfBranches()->Remove(qz);
    T->GetListOfBranches()->Remove(qv);

    TLeaf *lqx = T->GetLeaf("Q_Reco_X");
    TLeaf *lqy = T->GetLeaf("Q_Reco_Y");
    TLeaf *lqz = T->GetLeaf("Q_Reco_Z");
    TLeaf *lqv = T->GetLeaf("Q_Fit_Valid");
    T->GetListOfLeaves()->Remove(lqx);
    T->GetListOfLeaves()->Remove(lqy);
    T->GetListOfLeaves()->Remove(lqz);
    T->GetListOfLeaves()->Remove(lqv);

    T->Write();
    f->Write();
    f->Close();
  }
}
