/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.30
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCCOMPONENT_H
#define DCCOMPONENT_H
#include "BasicTypes.h"

//Class for Deferred Coupling CBMC components

namespace cbmc
{
class TrialMol;

class DCComponent
{
public:
  //Perform Decoupled portions of CBMC
  virtual void PrepareNew(TrialMol& newMol, uint molIndex) {}
  virtual void PrepareOld(TrialMol& oldMol, uint molIndex) {}

  //Perform Coupled final build
  virtual void BuildOld(TrialMol& oldMol, uint molIndex) = 0;
  virtual void BuildNew(TrialMol& newMol, uint molIndex) = 0;

  virtual void UpdateAcceptance(const TrialMol& mol) {}
  virtual ~DCComponent() {};
};
}

#endif
