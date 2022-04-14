#ifndef __MergerChainFactory__
#define __MergerChainFactory__

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

class MergerTChain;
class MergerTFile;

class MergerChainFactory
{
  public:
    MergerChainFactory( MergerConfig* _config, TRandom3* rndm, bool verbose );
    ~MergerChainFactory();

    // Methods
    std::vector<MergerTChain*> chainList;
    TRandom3* rndm;
    MergerConfig* config;
    MergerTChain* nextChain;
    std::map<double, int> timeComponentMap;
    std::pair<double, int> bufferTC;
    double bufferTimePrev;
    double bufferXPrev;
    double bufferYPrev;
    double bufferZPrev;
    int bufferFileIndex;
    int bufferEvtIndex;
    // cuts
    double time_window;
    double pos_window;

    double doubleBufferXPrev;
    double doubleBufferYPrev;
    double doubleBufferZPrev;
    double timenow;
    bool verbose;

    std::string LastFileName;

    // Member functions
    double nextEvent();
    void buildNewFile(std::string fname);
    std::vector<std::string> listDir(std::string directory);
};

class MergerTChain
{
  public:
    MergerTChain( std::string dstree, std::string dsbranch, std::string name, 
        std::vector<std::string> directory, double rate, bool is_single, TRandom3* rndm );
    ~MergerTChain();

    std::string dstree;
    std::string dsbranch;
    std::string name;
    std::vector<std::string> directory;
    bool is_single;
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

    std::vector<MergerTFile*> dataVec;
    std::vector<std::vector<double>> xpos;
    std::vector<std::vector<double>> ypos;
    std::vector<std::vector<double>> zpos;

    void addTime(double, int, int);
    void eventBuilder(bool);
    void shuffleDS();
    void addNewFile( std::string fname );
    void reset();
    void setupHeader();
    void setupDB();
    void getRandomEvent();
};

class MergerTFile
{
  public:
    MergerTFile( std::string dstree, std::string dsbranch, std::string fname, TRandom3* rndm );
    ~MergerTFile();

    std::string dstree;
    std::string dsbranch;
    std::string fname;
    TRandom3* rndm;

    TFile* tfile;
    TTree* ttree;
    RAT::DS::Root* ds;
    int entries;

    std::vector<RAT::DS::Root> getSubset(std::vector<int>);
    bool checkEvent();
    void open();
    void close();

};

#endif
