#include "FFDihedrals.h"  //Parent class
#include "FFSetup.h" //For initialization data
#include <algorithm> //for vector copying


void FFDihedrals::Init(ff_setup::Dihedral const& dih)
{
  uint size = dih.getTerms(), numSubDiv = dih.getnamecnt(), count = 0;
  Kchi = new double[size];
  n = new uint[size];
  delta = new double[size];
  subdiv.Init(numSubDiv);
  for (uint s = 0; s < numSubDiv; s++)
  {
    std::string div = dih.getname(s);
    uint cnt = dih.append(div, Kchi, delta, n, count);
    subdiv.Set(s, count, cnt);
    count += cnt;
  }
}

FFDihedrals::~FFDihedrals()
{
  delete[] Kchi;
  delete[] delta;
  delete[] n;
}