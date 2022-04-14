#include <MicroDS.hh>
#include <TTree.h>

MicroDS::MicroDS(){
  pdgcodes = new std::vector<Int_t>;
  mcKEnergies = new std::vector<double>;
  mcPosx = new std::vector<double>;
  mcPosy = new std::vector<double>;
  mcPosz = new std::vector<double>;
  mcDirx = new std::vector<double>;
  mcDiry = new std::vector<double>;
  mcDirz = new std::vector<double>;

  // Reconstructed variables / ev
  pedestal = new std::vector<int>;
  n100     = new std::vector<int>;
  n400     = new std::vector<int>;
  nReal    = new std::vector<int>;
  // Position
  x = new std::vector<double>;
  y = new std::vector<double>;
  z = new std::vector<double>;
  // Direction
  u = new std::vector<double>;
  v = new std::vector<double>;
  w = new std::vector<double>;
  // EV Info
  subev = new std::vector<int>;
  triggertime = new std::vector<double>;

}

MicroDS::~MicroDS()
{
  cleardata();
  delete pdgcodes;
  delete mcKEnergies;
  delete mcPosx;
  delete mcPosy;
  delete mcPosz;
  delete mcDirx;
  delete mcDiry;
  delete mcDirz;
  delete pedestal;
  delete n100;
  delete n400;
  delete nReal;
  delete x;
  delete y;
  delete z;
  delete u;
  delete v;
  delete w;
  delete subev;
  delete triggertime;
}

MicroDS MicroDS::clone()
{
  MicroDS mds;
  mds.mcT = mcT;
  mds.mcpcount = mcpcount;
  for(auto iv : *pdgcodes) mds.pdgcodes->push_back( iv );
  for(auto iv : *mcKEnergies) mds.mcKEnergies->push_back( iv );
  for(auto iv : *mcPosx) mds.mcPosx->push_back( iv );
  for(auto iv : *mcPosy) mds.mcPosy->push_back( iv );
  for(auto iv : *mcPosz) mds.mcPosz->push_back( iv );
  for(auto iv : *mcDirx) mds.mcDirx->push_back( iv );
  for(auto iv : *mcDiry) mds.mcDiry->push_back( iv );
  for(auto iv : *mcDirz) mds.mcDirz->push_back( iv );
  mds.evcount = evcount;
  for(auto iv : *pedestal) mds.pedestal->push_back( iv );
  for(auto iv : *n100) mds.n100->push_back( iv );
  for(auto iv : *n400) mds.n400->push_back( iv );
  for(auto iv : *nReal) mds.nReal->push_back( iv );
  for(auto iv : *x) mds.x->push_back( iv );
  for(auto iv : *y) mds.y->push_back( iv );
  for(auto iv : *z) mds.z->push_back( iv );
  for(auto iv : *u) mds.u->push_back( iv );
  for(auto iv : *v) mds.v->push_back( iv );
  for(auto iv : *w) mds.w->push_back( iv );
  for(auto iv : *subev) mds.subev->push_back( iv );
  for(auto iv : *triggertime) mds.triggertime->push_back( iv );
  return mds;
}

void MicroDS::cleardata()
{
    // MCParticles
    mcpcount = 0;
    pdgcodes->clear();
    mcKEnergies->clear();
    mcPosx->clear();
    mcPosy->clear();
    mcPosz->clear();
    mcDirx->clear();
    mcDiry->clear();
    mcDirz->clear();

    // Reconstructed variables / ev
    evcount;
    pedestal->clear();   // (-150, -50)
    n100->clear();       // (-20, 80)
    n400->clear();       // (-50, 350)
    nReal->clear();      // (-200, 600) Excluding noise hits
    // Position
    x->clear();
    y->clear();
    z->clear();
    // Direction
    u->clear();
    v->clear();
    w->clear();
    // EV Info
    subev->clear();
    triggertime->clear();
}

void MicroDS::NewBranches( TTree* t )
{
  t->Branch("mcpcount", &mcpcount);
  // Vector branches
  t->Branch("pdg", &pdgcodes);
  t->Branch("mcKEnergy", &mcKEnergies);
  t->Branch("mcposx", &mcPosx);
  t->Branch("mcposy", &mcPosy);
  t->Branch("mcposz", &mcPosz);
  t->Branch("mcposx", &mcDirx);
  t->Branch("mcposy", &mcDiry);
  t->Branch("mcposz", &mcDirz);
  // Vector EV
  t->Branch("evcount", &evcount);
  t->Branch("subev", &subev);
  t->Branch("triggertime", &triggertime);
  t->Branch("pedestal", &pedestal);
  t->Branch("n100", &n100);
  t->Branch("n400", &n400);
  t->Branch("nReal", &nReal);
  t->Branch("x", &x);
  t->Branch("y", &y);
  t->Branch("z", &z);
  t->Branch("u", &u);
  t->Branch("v", &v);
  t->Branch("w", &w);
}

void MicroDS::SetBranches( TTree* t )
{
  t->SetBranchAddress("mcpcount", &mcpcount);
  // Vector branches
  t->SetBranchAddress("pdg", &pdgcodes);
  t->SetBranchAddress("mcKEnergy", &mcKEnergies);
  t->SetBranchAddress("mcposx", &mcPosx);
  t->SetBranchAddress("mcposy", &mcPosy);
  t->SetBranchAddress("mcposz", &mcPosz);
  t->SetBranchAddress("mcposx", &mcDirx);
  t->SetBranchAddress("mcposy", &mcDiry);
  t->SetBranchAddress("mcposz", &mcDirz);
  // Vector EV
  t->SetBranchAddress("evcount", &evcount);
  t->SetBranchAddress("subev", &subev);
  t->SetBranchAddress("triggertime", &triggertime);
  t->SetBranchAddress("pedestal", &pedestal);
  t->SetBranchAddress("n100", &n100);
  t->SetBranchAddress("n400", &n400);
  t->SetBranchAddress("nReal", &nReal);
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);
  t->SetBranchAddress("u", &u);
  t->SetBranchAddress("v", &v);
  t->SetBranchAddress("w", &w);
}
