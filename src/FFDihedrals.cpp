#include "FFDihedrals.h"  //Parent class
#include "FFSetup.h" //For initialization data
#include <algorithm> //for vector copying


void FFDihedrals::Init(ff_setup::Dihedral const& dih)
{
   uint size = dih.countTerms, numSubDiv = dih.name.size(), count=0;
   Kchi = new double [size];
   n = new uint [size];
   delta = new double [size];
   subdiv.Init(numSubDiv);
   for (uint s = 0; s < numSubDiv; s++)
   {
      std::string div = dih.name[s];
      std::map<std::string, std::vector<uint> >::const_iterator itUInt =
         dih.n.find(div);
      std::map<std::string, std::vector<double> >::const_iterator itDbl = 
         dih.Kchi.find(div);
      subdiv.Set(s, count, itUInt->second.size());
      std::copy(itDbl->second.begin(), itDbl->second.end(), Kchi+count);
      std::copy(itUInt->second.begin(), itUInt->second.end(), n+count);
      itDbl = dih.delta.find(div);
      std::copy(itDbl->second.begin(), itDbl->second.end(), delta+count);
      count += itUInt->second.size();
   }
}

FFDihedrals::~FFDihedrals()
{ 
   delete[] Kchi; 
   delete[] delta; 
   delete[] n;
}

