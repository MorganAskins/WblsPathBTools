#include <MergerChainFactory.hh>
#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <utility>
#include <boost/filesystem.hpp>
#include <TChain.h>
#include <TTimeStamp.h>

// Factory manages each component of the model, pointing to a list of Chains
// Each Chain then handles the list of files
//
// Key points:
// 1. Root files love to be read sequentially
// 2. Root files love to be opened and read in order

MergerChainFactory::MergerChainFactory( MergerConfig* _config, TRandom3* rndm, bool verbose ) :
  config(_config), rndm(rndm), verbose(verbose)
{
  for( auto mcc : config->componentList )
  {
    // Construct TChain and store in factory
    std::string chaindir = config->baseDir + "/" + mcc->dir;
    std::vector<std::string> rootfiles = listDir( chaindir );
    chainList.push_back( new MergerTChain( config->dstree, config->dsbranch, 
          mcc->name, rootfiles, mcc->rate, mcc->is_single, rndm ) );
    this->LastFileName = rootfiles[0];
  }
  this->bufferTC = std::make_pair(-100, 0);
  this->bufferTimePrev = -1000;
  this->bufferXPrev = 9999999;
  this->bufferYPrev = 9999999;
  this->bufferZPrev = 9999999;
  this->doubleBufferXPrev = 2*bufferXPrev;
  this->doubleBufferYPrev = 2*bufferYPrev;
  this->doubleBufferZPrev = 2*bufferZPrev;
  this->bufferFileIndex = 0;
  this->bufferEvtIndex = 0;
  this->timenow = 0;
  this->time_window = config->deltat;
  this->pos_window = config->deltar;
  for( auto cl : chainList )
  {
    printf("%f %f\n", cl->rate, cl->efficiency);
  }
}

MergerChainFactory::~MergerChainFactory()
{
  // Delete containers
  for(auto p : this->chainList) delete p;
  this->chainList.clear();
}

double MergerChainFactory::nextEvent()
{
  // Choose what the next event will be, but do not access it yet
  // We must build a map of < component, vector<time> >
  std::map<double, MergerTChain*> order;
  std::map<double, int> timeIndex;
  // Loop through the chain and get the "next" event of each.
  for(int cl=0; cl<chainList.size(); cl++)
  {
    auto mtc = chainList[cl];
    double u = this->rndm->Rndm();
    double poissontime     = -log(1-u)/( mtc->rate * mtc->efficiency );
    order[poissontime]     = mtc;
    timeIndex[poissontime] = cl;
  }
  // Of the "next" possible events, choose the soonest
  nextChain = order.begin()->second;
  timenow += order.begin()->first;
  // Decide what to do with the last buffer (write or throw)
  // << Timing information
  double lookback    = this->bufferTC.first - bufferTimePrev;
  double lookforward = timenow - this->bufferTC.first;
  // << Position information
  nextChain->getRandomEvent(); // This sets the MergerTChain file_index and evt_index
  double pforward = sqrt( pow( nextChain->x - bufferXPrev, 2 ) +
                          pow( nextChain->y - bufferYPrev, 2 ) +
                          pow( nextChain->z - bufferZPrev, 2 ) );
  double pback    = sqrt( pow( bufferXPrev - doubleBufferXPrev, 2 ) +
                          pow( bufferYPrev - doubleBufferYPrev, 2 ) +
                          pow( bufferZPrev - doubleBufferZPrev, 2 ) );
  // Check the last event for writing
  // if ( 
  //     (
  //      ( (lookback < time_window) || (lookforward < time_window) ) && //Time Cut
  //      ( (pforward < pos_window) || (pback < pos_window) ) //Position cut
  //     ) || ( !chainList[this->bufferTC.second]->is_single ) // We should also ignore this if the event is a multi
  //    )
  if ( ( ( (lookback < time_window) && (pback < pos_window) ) ||
         ( (lookforward < time_window) && (pforward < pos_window) ) ) ||
       ( !chainList[this->bufferTC.second]->is_single ) )
  {
    chainList[ this->bufferTC.second ]->addTime( this->bufferTC.first, 
        bufferFileIndex, bufferEvtIndex );
    // printf("Pos: %f, %f -- DT: %f, %f\n", 
    //     pforward, pback, lookback, lookforward);
    // printf("\tFile: %i, Evt: %i\n", bufferFileIndex, bufferEvtIndex);
    // printf("\t(%f, %f, %f)\n", bufferXPrev, bufferYPrev, bufferZPrev);
    // New method should be addEvent( time, file_index, evt_index );
    timeComponentMap.insert( bufferTC );
  }
  // Write this event to the buffer
  bufferTimePrev = this->bufferTC.first;
  this->bufferTC = std::make_pair( timenow, timeIndex.begin()->second );
  // Update two events ago
  doubleBufferXPrev = bufferXPrev;
  doubleBufferYPrev = bufferYPrev;
  doubleBufferZPrev = bufferZPrev;
  // Update Buffer
  bufferXPrev     = nextChain->x;
  bufferYPrev     = nextChain->y;
  bufferZPrev     = nextChain->z;
  bufferFileIndex = nextChain->file_index;
  bufferEvtIndex  = nextChain->evt_index;
  return timenow;
}

void MergerChainFactory::buildNewFile(std::string fname)
{
  // Build vectors of ds events on each chain
  for( auto mtc : chainList )
  {
    mtc->eventBuilder(verbose=this->verbose);
  }
  // Top file
  std::cout << LastFileName << std::endl;
  TFile* oldFile = new TFile(LastFileName.c_str());
  // TFile* oldFile = nextChain->dataVec[0]->tfile;
  TTree* oldRunTree = (TTree*)oldFile->Get("runT");

  // Write to file
  std::string outname = config->trainingDir + "/" + fname;
  TFile* f = new TFile(outname.c_str(), "recreate");
  // Lets add a special header with info from this merge
  TTree* header = new TTree("header", "Merger information");
  double time = timenow; // seconds I believe
  header->Branch("livetime", &time);
  header->Fill();
  // Run tree as well, we could clone from a random file?
  TTree* runT = oldRunTree->CloneTree(0);
  oldRunTree->GetEvent(0);
  runT->Fill();

  TTree* t = new TTree("T", "merged");
  RAT::DS::Root* ds = new RAT::DS::Root();
  t->Branch("ds", &ds);
  std::string iname;
  t->Branch("name", &iname);
  // Combine the chains into a single file
  if(verbose)
    printf("Writing to file %s ...", outname.c_str());
  for( auto iv : this->timeComponentMap )
  {
    double time = iv.first;
    int idx = iv.second;
    MergerTChain* mtc = chainList[idx];
    iname = mtc->name;
    ds = &( *(mtc->dsitr) );
    ++mtc->dsitr;
    if( ds->ExistMC() )
    {
      // Update the event. Set simulation time to Jan 1st 1970.
      // Beware ... root sucks ...
      // Also, uses 32 bit int, so TTimeStamp dies in 1938 -.-
      TTimeStamp tt(1970, 1, 1, 0, 0, 0);
      // Grab the MC time
      time_t seconds = static_cast<time_t>(floor(time));
      Int_t nanoseconds = static_cast<Int_t>( (time - seconds)*1e9 );
      TTimeStamp mctime(seconds, nanoseconds);
      //std::cout << "SET: " << mctime.GetSec() << " & " << nanoseconds << std::endl;
      //mctime.Add(tt);
      // 
      RAT::DS::MC* mc = ds->GetMC();
      mc->SetUTC( mctime );
      // Print some things here, delete after
      TVector3 v = ds->GetEV(0)->GetPathFit()->GetPosition();
      //printf("True pos: %f, %f, %f\n", v.X(), v.Y(), v.Z());
      //
      t->Fill();
    }
  }
  if( verbose )
    printf(" done\n");
  f->Write(0, TObject::kOverwrite);
  f->Close();

  delete f;

  for( auto mtc : chainList )
  {
    mtc->reset();
  }
  timeComponentMap.clear();
  this->timenow = 0;
}

std::vector<std::string> MergerChainFactory::listDir(std::string directory)
{
  std::vector<std::string> files;
  for(auto &p : boost::filesystem::directory_iterator( directory ))
  {
    files.push_back( p.path().string() );
  }
  return files;
}

// Merger TChain
MergerTChain::MergerTChain( std::string dstree, std::string dsbranch, 
    std::string name, std::vector<std::string> directory, double rate, 
    bool is_single, TRandom3* rndm ) :
  dstree(dstree), dsbranch(dsbranch), name(name), directory(directory), 
  rate(rate), is_single(is_single), rndm(rndm)
{
  ds = new RAT::DS::Root();
  setupHeader();
  setupDB();
}


void MergerTChain::setupHeader()
{
  // Look in the header chain for the efficiency
  this->counter = 0;
  this->efficiency = 0.0;
  double subEfficiency;
  int count=0;
  int totalcount = directory.size();
  for(auto v : directory )
  {
    addNewFile( v );
    TFile* f = TFile::Open(v.c_str());
    TTree* headchain = (TTree*)f->Get("header");
    headchain->SetBranchAddress("efficiency", &subEfficiency);
    headchain->GetEvent(0);
    if( isnan(subEfficiency) )
      subEfficiency = 0.0;
    this->efficiency += subEfficiency / totalcount;
    f->Close();
    delete f;
    printf("Header -> %i / %i -- (%f)\r", count, totalcount, subEfficiency);
    count++;
  }
  printf("\n\teff: %f\n", this->efficiency );
  this->entries = dataVec.size();
  printf("\nEntries: %i\n", this->entries);
}

void MergerTChain::setupDB()
{
  // Database TChain
  int count = 0;
  int totalcount = directory.size();
  for(auto v : directory )
  {
    TFile* f = TFile::Open(v.c_str());
    TTree* dbchain = (TTree*)f->Get("posdb");
    std::vector<double>* vxpos = new std::vector<double>;
    std::vector<double>* vypos = new std::vector<double>;
    std::vector<double>* vzpos = new std::vector<double>;
    dbchain->SetBranchAddress("xdb", &vxpos);
    dbchain->SetBranchAddress("ydb", &vypos);
    dbchain->SetBranchAddress("zdb", &vzpos);
    dbchain->GetEvent(0);
    this->xpos.push_back( *vxpos );
    this->ypos.push_back( *vypos );
    this->zpos.push_back( *vzpos );
    f->Close();
    delete f;
    delete vxpos;
    delete vypos;
    delete vzpos;
    printf("DB-> %i / %i\r", count, totalcount);
    count++;
  }
  printf("\n");
}

void MergerTChain::getRandomEvent()
{
  // Choose a random file
  file_index = int( rndm->Rndm() * entries );
  // Choose a random evt from the file
  evt_index  = int( rndm->Rndm() * xpos[ file_index ].size() );
  x = xpos[file_index][evt_index];
  y = ypos[file_index][evt_index];
  z = zpos[file_index][evt_index];
}

void MergerTChain::addNewFile( std::string fname )
{
  MergerTFile* mtf = new MergerTFile( dstree, dsbranch, fname, rndm );
  this->dataVec.push_back(mtf);
}

MergerTChain::~MergerTChain()
{
  // Lots to do here
  for(auto p : this->dataVec) delete p;
  this->dataVec.clear();
}

void MergerTChain::eventBuilder(bool verbose=false)
{
  // DEBUG
  // for( auto v : evtStamps )
  //   printf("%i _ ", v);
  // printf("\n");
  // Info stored in timeStamps, fileStamps, evtStamps
  std::vector< std::vector<int> > fcount(entries);
  for(int i=0; i < fileStamps.size(); i++)
  {
    fcount[ fileStamps[i] ].push_back( evtStamps[i] );
  }
  // Sort
  for(auto &v : fcount)
  {
    std::sort( v.begin(), v.end() );
    // for(auto k : v)
    //   printf("%i ", k);
    // printf("\n");
  }

  std::vector<RAT::DS::Root> dsholder;

  for(int iv=0; iv < entries; iv++)
  {
    if( verbose )
      printf("\t<eventbuilder>: Files %i of %i\t\t(%i)\r", iv, entries, fcount[iv].size());
    if( fcount[iv].size() > 0 )
    {
      dataVec[iv]->open();
      std::vector<RAT::DS::Root> a = dataVec[iv]->getSubset( fcount[iv] );
      dsholder.insert( dsholder.end(), a.begin(), a.end() );
      dataVec[iv]->close();
      delete dataVec[iv]; // Does this help?
    }
  }
  if( verbose )
    printf("\n");
  // Shuffle vector<ds>
  //this->shuffleDS();
  dsevents.resize( fileStamps.size() );
  for( int i=0; i < fileStamps.size(); i++ )
  {
    int findex  = fileStamps[i];
    int evindex = evtStamps[i];
    // Where in dsholder do I live?
    auto v = fcount[findex];
    int index = std::distance( v.begin(), std::find(v.begin(), v.end(), evindex) );
    // printf("evtstamp: %i at %i\n", evindex, index);
    dsevents[i] = dsholder[index];
  }
  // fcount has the new [file_index][evt] order
  // Initialize / reset the event iterator
  this->dsitr = this->dsevents.begin();
}

void MergerTChain::shuffleDS()
{
  auto first = dsevents.begin();
  auto last = dsevents.end();
  for(auto i=(last-first)-1; i>0; --i)
  {
    int u = int( rndm->Rndm() * i );
    std::swap(first[i], first[u]);
  }
}

void MergerTChain::addTime(double t, int findex, int evindex)
{
  timeStamps.push_back(t);
  fileStamps.push_back(findex);
  evtStamps.push_back(evindex);
}

void MergerTChain::reset()
{
  timeStamps.clear();
  fileStamps.clear();
  evtStamps.clear();
}

// Control individual TFiles (lowest level)

MergerTFile::MergerTFile( std::string dstree, std::string dsbranch, std::string fname, TRandom3* rndm ) :
  dstree(dstree), dsbranch(dsbranch), fname(fname), rndm(rndm)
{
}

MergerTFile::~MergerTFile()
{
  //tfile->Close();
}

void MergerTFile::open()
{
  tfile = TFile::Open( fname.c_str(), "read" );
  ttree = (TTree*)tfile->Get("T");
  ds = new RAT::DS::Root();
  ttree->SetBranchAddress( "ds", &ds );
  entries = ttree->GetEntries();
}

void MergerTFile::close()
{
  tfile->Close();
  delete tfile;
  delete ds;
}

std::vector<RAT::DS::Root> MergerTFile::getSubset(std::vector<int> events)
{
  std::vector<RAT::DS::Root> ratpile;
  if( events.size() < 1 ) return ratpile;

  //std::vector<int> entry;
  //for(int i=0; i<num; i++)
  //{
  //  entry.push_back( int( rndm->Rndm()*entries ) );
  //}
  std::sort(events.begin(), events.end());
  for(auto iv : events)
  {
    ttree->GetEvent(iv);
    if( this->checkEvent() )
    {
      ratpile.push_back(*ds);
    }
    else
    {
      ratpile.push_back(RAT::DS::Root());
    }
    //ratpile.push_back(*ds);
  }
  // Did it work??
  return ratpile;
}

bool MergerTFile::checkEvent()
{
  // Not using this yet
  return true;
  // Right now this is just an example, will change based on Stephane's structure.
  int min_nhit = 60;
  // High level cuts on event ev branch
  // 1. EVCount < 1
  if( ds->GetEVCount() < 1 )
    return false;
  // 2. nhits < min_nhit
  int high_count = 0;
  for( int iev=0; iev < ds->GetEVCount(); iev++ )
  {
    RAT::DS::EV* ev = ds->GetEV(iev);
    int count = ev->GetPMTCount();
    if( count > high_count )
      high_count = count;
  }
  if( high_count < min_nhit )
    return false;

  return true;
}
