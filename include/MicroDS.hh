#ifndef __MicroDS__
#define __MicroDS__

#include <TTimeStamp.h>
#include <TTree.h>
#include <vector>

class MicroDS
{
  public:
    MicroDS();
    ~MicroDS();

    void cleardata();
    void NewBranches(TTree*);
    void SetBranches(TTree*);
    MicroDS clone();

    TTimeStamp mcT;
    // MCParticles
    int mcpcount;
    std::vector<Int_t>* pdgcodes;
    std::vector<double>* mcKEnergies;
    std::vector<double>* mcPosx;
    std::vector<double>* mcPosy;
    std::vector<double>* mcPosz;
    std::vector<double>* mcDirx;
    std::vector<double>* mcDiry;
    std::vector<double>* mcDirz;

    // Reconstructed variables / ev
    int evcount;
    std::vector<int>* pedestal;   // (-150, -50)
    std::vector<int>* n100;       // (-20, 80)
    std::vector<int>* n400;       // (-50, 350)
    std::vector<int>* nReal;      // (-200, 600) Excluding noise hits
    // Position
    std::vector<double>* x;
    std::vector<double>* y;
    std::vector<double>* z;
    // Direction
    std::vector<double>* u;
    std::vector<double>* v;
    std::vector<double>* w;
    // EV Info
    std::vector<int>*        subev;
    std::vector<double>* triggertime;

};

#endif
