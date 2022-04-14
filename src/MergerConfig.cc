#include <MergerConfig.hh>
// Boost-json
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>

MergerConfig::MergerConfig( std::string config_file, std::string subdir="" ) :
  configFile(config_file)
{
  boostRead();
  if( subdir != "" )
  {
    this->baseDir += "/" + subdir;
    this->trainingDir += "/" + subdir;
    this->sampleDir += "/" + subdir;
  }
}

MergerConfig::~MergerConfig()
{
  // Delete containers
  for(auto p : this->componentList) delete p;
  this->componentList.clear();
}

void MergerConfig::boostRead()
{
  // property tree from file
  namespace pt = boost::property_tree;
  pt::ptree iroot;

  // Load the json file
  pt::read_json( this->configFile, iroot );
  // Grab the top level dictionaries
  std::string component_list = "components";
  
  // Header information
  pt::ptree header  = iroot.get_child( "header" );
  this->baseDir     = header.get<std::string>( "base_directory" );
  this->trainingDir = header.get<std::string>( "training_directory" );
  this->sampleDir   = header.get<std::string>( "sample_directory" );
  this->dstree      = header.get<std::string>( "dstree" );
  this->dsbranch    = header.get<std::string>( "dsbranch" );
  this->deltat      = header.get<double>( "deltat" );
  this->deltar      = header.get<double>( "deltar" );

  // Grab each component and store its name, directory, and rate
  for( const auto& parent : iroot.get_child(component_list) )
  {
    pt::ptree component = parent.second;
    std::string pname = component.get<std::string>("name");
    this->componentList.push_back( new MCComponent( 
          component.get<std::string>("name"),
          component.get<std::string>("directory"),
          component.get<double>("rate"),
          component.get<std::string>("class")
          ));
  }
}

void MergerConfig::print()
{
  std::cout << "MergerConfig (" << this->configFile << ")" << std::endl;
  std::cout << "\t" << "Base Directory     :" << this->baseDir << std::endl;
  std::cout << "\t" << "Training Directory :" << this->trainingDir << std::endl;
  std::cout << "\t" << "Sample Directory   :" << this->sampleDir << std::endl;
  std::cout << "\t" << "Read Tree          :" << this->dstree << std::endl;
  std::cout << "\t" << "Read Branch        :" << this->dsbranch << std::endl;
  for( auto mcc : componentList )
  {
    std::cout << "\t" << "> " << mcc->name << " @ " << mcc->rate << " :" << mcc->dir << std::endl;
  }
}
