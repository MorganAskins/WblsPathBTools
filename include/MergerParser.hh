#ifndef __MergerParser__
#define __MergerParser__

#include <string>
#include <vector>

class MergerParser
{
  public:
    MergerParser( std::vector<std::string> _args);
    ~MergerParser();

    std::vector<std::string> args;

    // Arguments
    std::string config;
    int num;
    int start;
    double time;
    bool verbose;
    bool superverbose;
    std::string subdir;
  private:
    void help();
    void setDefaultParams();
};

#endif
