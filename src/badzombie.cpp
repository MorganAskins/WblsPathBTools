#include <TFile.h>
#include <cstdio>
#include <string>
#include <vector>
using namespace std;

void destroy_all_zombies( vector<string> zombies );

int main(int argc, char** argv)
{
  vector<string> files(argv+1, argv+argc);
  vector<string> zombies;
  for( auto f : files )
  {
    TFile xf(f.c_str(), "Q");
    if( xf.IsZombie() )
    {
      zombies.push_back(f);
    }
    else if( xf.GetNkeys() < 1 )
    {
      zombies.push_back(f);
    }
    else if( ! xf.GetListOfKeys()->Contains("T") )
    {
      zombies.push_back(f);
    }
    xf.Close();
  }
  destroy_all_zombies( zombies );
  printf("Shaun has smashed %i zombies with a cricket bat\n", zombies.size());
  return 0;
}

void destroy_all_zombies( vector<string> zombies )
{
  for( auto f : zombies )
  {
    remove(f.c_str());
  }
}
