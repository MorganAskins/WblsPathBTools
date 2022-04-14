#ifndef __MergerConfig__
#define __MergerConfig__

#include <string>
#include <vector>

class MCComponent;

class MergerConfig
{
  public:
    MergerConfig( std::string config_file, std::string subdir );
    ~MergerConfig();
    void boostRead();
    void print();

    std::vector<MCComponent*> componentList;
    std::string configFile;
    std::string baseDir;
    std::string trainingDir;
    std::string sampleDir;
    std::string dstree;
    std::string dsbranch;
    double deltat;
    double deltar;
};

class MCComponent
{
  public:
    MCComponent( std::string _name, std::string _dir, double _rate, std::string classify ) :
      name(_name), dir(_dir), rate(_rate) {
        if( classify == "single" )
          is_single = true;
        else
          is_single = false;
      };

    std::string name;
    std::string dir;
    double rate;
    bool is_single;
};

#endif
