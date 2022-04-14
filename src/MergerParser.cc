#include <MergerParser.hh>
// Boost-json
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <stdlib.h>

MergerParser::MergerParser( std::vector<std::string> _args) :
  args(_args)
{
  // Help > all other commands and exits
  if( args.size() == 0 ) help();
  for( auto v : args )
  {
    if( v == "--help" || v == "-h" ) help();
  }
  setDefaultParams();

  // Json config file
  this->config = args[0];
  std::string iv = "";
  for( auto v : args )
  {
    // Number of datasets to produce
    if( iv == "-n" || iv == "--num" )
    {
      this->num = stoi(v);
    }
    // Where to start counting from (useful in batch mode)
    if( iv == "-s" || iv == "--start" )
    {
      this->start = stoi(v);
    }
    // Time in seconds of the dataset
    if( iv == "-t" || iv == "--time" )
    {
      this->time = stod(v);
    }
    // Subdirectory
    if( iv == "-d" || iv == "--subdir" )
    {
      this->subdir = v;
    }
    // Verbose
    if( v == "-v" || v == "--verbose" )
    {
      this->verbose = true;
    }
    // REALLY verbose
    if( v == "-vv" )
    {
      this->superverbose = true;
    }
    iv = v;
  }
  // Verbose print
  if( this->verbose )
  {
    std::cout << "| ---------------------------------------------" << std::endl;
    std::cout << "| MergerParser" << std:: endl;
    std::cout << "| Configuration  : " << this->config << std::endl;
    std::cout << "| Num datasets   : " << this->num << std::endl;
    std::cout << "| Dataset start  : " << this->start << std::endl;
    std::cout << "| Dataset length : " << this->time << std::endl;
    std::cout << "| Subdirectory   : " << this->subdir << std::endl;
    std::cout << "| ---------------------------------------------" << std::endl;
  }
}

MergerParser::~MergerParser()
{
}

void MergerParser::setDefaultParams()
{
  this->num     = 1;
  this->time    = 3600;
  this->verbose = false;
  this->subdir  = "wm_20pct_geo/wbls_1pct";
}

void MergerParser::help()
{
  std::cout << "mergeddatasets <configfile.json> <options>" << std::endl;
  std::cout << "    -h,--help    : Print help dialog" << std::endl;
  std::cout << "    -n,--num     : Specify number of datasets" << std::endl;
  std::cout << "    -t,--time    : Length of dataset (seconds)" << std::endl;
  std::cout << "    -d,--subdir  : Subdirectory (geo/target)" << std::endl;
  std::cout << "    -v,--verbose : Verbose" << std::endl;
  exit(EXIT_SUCCESS);
}
