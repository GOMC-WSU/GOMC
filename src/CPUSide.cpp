/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "System.h"
#include "StaticVals.h"
#include "CPUSide.h" //Spec declaration
#include "GOMCEventsProfile.h"


CPUSide::CPUSide(System & sys, StaticVals & statV, Setup & set) :
  varRef(sys, statV, set.mol.molVars.moleculeKindNames), pdb(sys, statV), 
  psf(statV.mol, sys, set), 
  xstBinary(sys, statV),
  console(varRef), block(varRef),
  hist(varRef), checkpoint(sys, statV)
#if ENSEMBLE == GCMC
  , sample_N_E(varRef)
#endif
#if ENSEMBLE == NVT || ENSEMBLE == NPT
  , freeEnergy(varRef, sys)
#endif
{}

void CPUSide::Init( PDBSetup const& pdbSet, 
                    config_setup::Output const& out,
                    const ulong tillEquil, 
                    const ulong totSteps, 
                    ulong startStep)
{
  equilSteps = tillEquil;
  //Initialize arrays in object that collects references and calc'ed vals.
  varRef.Init();
  //Initialize output components.
  timer.Init(out.console.frequency, totSteps, startStep);
  outObj.push_back(&console);
  outObj.push_back(&pdb);
  outObj.push_back(&xstBinary);
  outObj.push_back(&psf);
  if (out.statistics.settings.block.enable)
    outObj.push_back(&block);
  if (out.restart.settings.enable)
    outObj.push_back(&checkpoint);

#if ENSEMBLE == GCMC
  outObj.push_back(&hist);
  outObj.push_back(&sample_N_E);
#endif
#if ENSEMBLE == NVT || ENSEMBLE == NPT
  outObj.push_back(&freeEnergy);
#endif
  //Calculate pressure, heat of vap. (if applicable), etc.
  varRef.CalcAndConvert(0);
  for (uint o = 0; o < outObj.size(); o++)
    outObj[o]->Init(pdbSet.atoms, out, tillEquil, totSteps);
}

void CPUSide::Output(const ulong step)
{
  GOMC_EVENT_START(1, GomcProfileEvent::FILE_IO_OUTPUT);
  //Calculate pressure, heat of vap. (if applicable), etc.
  varRef.CalcAndConvert(step);
  //Do standard output events.
  for (uint o = 0; o < outObj.size(); o++)
    outObj[o]->Output(step);
  timer.CheckTime(step);
  GOMC_EVENT_STOP(1, GomcProfileEvent::FILE_IO_OUTPUT);
}
