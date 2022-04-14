#ifndef __MergerMicroChainFactory__
#define __MergerMicroChainFactory__

#include <MergerConfig.hh>
#include <MicroDS.hh>
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

class MergerMicroChain;
class MergerMicroFile;

class MergerMicroChainFactory
{
  public:
    MergerMicroChainFactory( MergerConfig* _config, TRandom3* rndm );
    ~MergerMicroChainFactory();

    // Methods
    std::vector<MergerMicroChain*> chainList;
    TRandom3* rndm;
    MergerConfig* config;
    MergerMicroChain* nextChain;
    std::map<double, int> timeComponentMap;
    double timenow;

    // Member functions
    double nextEvent();
    void buildNewFile(std::string fname);
    std::vector<std::string> listDir(std::string directory);
};

class MergerMicroChain
{
  public:
    MergerMicroChain( std::string dstree, std::string name, std::vector<std::string> directory, double rate, TRandom3* rndm );
    ~MergerMicroChain();

    std::string dstree;
    std::string name;
    std::vector<std::string> directory;
    double rate;
    double efficiency;
    int entries;
    int counter;

    std::vector<double> timeStamps;
    TRandom3* rndm;

    std::vector<MicroDS> dsevents;
    std::vector<MicroDS>::iterator dsitr;

    std::vector<MergerMicroFile*> dataVec;

    void addTime(double t);
    void eventBuilder();
    void shuffleDS();
    void addNewFile( std::string fname );
    void reset();
};

class MergerMicroFile
{
  public:
    MergerMicroFile( std::string dstree, std::string fname, TRandom3* rndm );
    ~MergerMicroFile();

    std::string dstree;
    std::string fname;
    TRandom3* rndm;

    TFile* tfile;
    TTree* ttree;
    MicroDS* ds;
    int entries;

    std::vector<MicroDS> getSubset(int num);
    bool checkEvent();
    void open();
    void close();
};

#endif
