#include "System.h"
#include "StaticVals.h"
#include "CPUSide.h" //Spec declaration

CPUSide::CPUSide(System & sys, StaticVals & statV) :
   varRef(sys, statV), pdb(sys, statV), console(varRef), block(varRef), 
   fluct(varRef), hist(varRef)
#if ENSEMBLE == GCMC
   , sample_N_E(varRef)
#endif
  {}

void CPUSide::Init(PDBSetup const& pdbSet, config_setup::Output const& out,
                   const ulong tillEquil, const ulong totSteps) 
{
   //Initialize arrays in object that collects references and calc'ed vals.
   varRef.Init(pdbSet.atoms);
   //Initialize output components.
   timer.Init(out.console.frequency, totSteps);
   outObj.push_back(&console);
   outObj.push_back(&pdb);
   outObj.push_back(&block);
   outObj.push_back(&fluct);
   outObj.push_back(&hist);
#if ENSEMBLE == GCMC
   outObj.push_back(&sample_N_E);
#endif
   //Calculate pressure, heat of vap. (if applicable), etc.
   varRef.CalcAndConvert();
   for (uint o = 0; o < outObj.size(); o++)
      outObj[o]->Init(pdbSet.atoms, out, tillEquil, totSteps);
}

void CPUSide::Output(const ulong step)
{
   //Calculate pressure, heat of vap. (if applicable), etc.
   varRef.CalcAndConvert();
   //Do standard output events.
   for (uint o = 0; o < outObj.size(); o++)
      outObj[o]->Output(step);
   timer.CheckTime(step);
}
