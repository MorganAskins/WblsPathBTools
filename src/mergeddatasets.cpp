#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <ctime>
#include <MergerConfig.hh>
#include <MergerParser.hh>
#include <MergerChainFactory.hh>

#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

int main(int argc, char** argv)
{
  // Parse commands
  std::vector<std::string> argument_vector(argv+1, argv+argc);
  MergerParser parser( argument_vector );

  // Read the config.json file to get event types, locations, and rates
  MergerConfig* config = new MergerConfig( parser.config, parser.subdir );
  if( parser.verbose ) config->print();

  TRandom3* rndm = new TRandom3();
  // Unique seeding time * pid
  long seed = time(nullptr) * getpid();
  rndm->SetSeed(seed);

  // Main loop
  for(int loop=parser.start; loop<(parser.num+parser.start); ++loop)
  {
    // Build input TChains
    MergerChainFactory factory( config, rndm, parser.superverbose );
    // Loop in time, grabbing entries based on poisson of rate
    double start_time = 0.0;
    if( parser.verbose )
      printf("Event: %i\n", loop);
    while( start_time < parser.time )
    {
      // if( parser.verbose )
      //   printf("\tTime: %f / %f\r", start_time, parser.time);
      double next_time = factory.nextEvent();
      start_time = next_time;
    }
    if( parser.verbose )
      printf("Total events: %i\n", factory.timeComponentMap.size());
    // File name
    std::stringstream ss;
    ss << "mergedfile_" << loop << ".root";
    std::string outfile_name = ss.str();
    if( parser.verbose )
      printf("::Writing out to %s\n", outfile_name.c_str());
    // Build data file
    factory.buildNewFile( outfile_name );
    std::cout << "Processed " << loop << " of " << parser.num << "\r" << std::flush;
  }

  delete rndm;
  delete config;
  return 0;
}
