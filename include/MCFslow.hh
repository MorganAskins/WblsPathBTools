#ifndef __MCFslow__
#define __MCFslow__

#include <MergerConfig.hh>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <string>
#include <vector>
#include <map>
#include <utility>

class MTCslow;
class MTFslow;

class MCFslow
{
  public:
    MCFslow( MergerConfig* _config, TRandom3* rndm, bool verbose );
    ~MCFslow();

    // Methods
    std::vector<MTCslow*> chainList;
    TRandom3* rndm;
    MergerConfig* config;
    MTCslow* nextChain;
    std::map<double, int> timeComponentMap;
    std::pair<double, int> bufferTC;
    double bufferTimePrev;
    double bufferXPrev;
    double bufferYPrev;
    double bufferZPrev;
    double bufferFileIndex;
    double bufferEvtIndex;

    double doubleBufferXPrev;
    double doubleBufferYPrev;
    double doubleBufferZPrev;
    double timenow;
    bool verbose;

    // Member functions
    double nextEvent();
    void buildNewFile(std::string fname);
    std::vector<std::string> listDir(std::string directory);
};

class MTCslow
{
  public:
    MTCslow( std::string dstree, std::string dsbranch, std::string name, std::vector<std::string> directory, double rate, TRandom3* rndm );
    ~MTCslow();

    std::string dstree;
    std::string dsbranch;
    std::string name;
    std::vector<std::string> directory;
    double rate;
    double efficiency;
    // Entries is the number of files
    int entries;
    int counter;
    // Current file and event index
    int file_index;
    int evt_index;
    double x, y, z;

    std::vector<double> timeStamps;
    std::vector<int> fileStamps;
    std::vector<int> evtStamps;
    TRandom3* rndm;

    RAT::DS::Root* ds;
    std::vector<RAT::DS::Root> dsevents;
    std::vector<RAT::DS::Root>::iterator dsitr;

    std::vector<MTFslow*> dataVec;
    std::vector<std::vector<double>> xpos;
    std::vector<std::vector<double>> ypos;
    std::vector<std::vector<double>> zpos;

    void addTime(double);
    void eventBuilder(bool);
    void shuffleDS();
    void addNewFile( std::string fname );
    void reset();
    void setupHeader();
    void setupDB();
    void getRandomEvent();
};

class MTFslow
{
  public:
    MTFslow( std::string dstree, std::string dsbranch, std::string fname, TRandom3* rndm );
    ~MTFslow();

    std::string dstree;
    std::string dsbranch;
    std::string fname;
    TRandom3* rndm;

    TFile* tfile;
    TTree* ttree;
    RAT::DS::Root* ds;
    int entries;

    std::vector<RAT::DS::Root> getSubset(int);
    bool checkEvent();
    void open();
    void close();

};

#endif
