#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TTimeStamp.h>
#include <TVector3.h>
#include <TH2F.h>

#include <RAT/DS/Root.hh>
#include <RAT/DSReader.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/PathFit.hh>
using namespace std;

#include <WbLS_Analyser_RATDS.C>

bool make_pdf(vector<string> p)
{
  for( auto v : p )
    if( v == "--pdf" )
      return true;
  return false;
}

string pdf_input(vector<string> p)
{
  // Grab the first argument that isn't --pdf
  for( auto v : p )
    if( v != "--pdf" )
      return v;
  return "NONE";
}

string new_name_and_file( string p )
{
  return p;
}

int main(int argc, char** argv)
{
  // Run as ./qfit --pdf input_files.root
  vector<string> args(argv+1, argv+argc);
  // Need to know if fitting or not
  if( make_pdf(args) )
  {
    // Generate the pdf files
    string pdf_file = pdf_input(args);
    if( pdf_file == "NONE" )
    {
      cout << "Error file not specified" << endl;
      exit(0);
    }
    cout << "Generating PDF" << endl;
    PDF_Producer( pdf_file.c_str() );
  }
  else
  {
    // Fit the given files
    for(auto v : args)
    {
      RATDS_Output( v.c_str(), 2 );
    }
  }
}
