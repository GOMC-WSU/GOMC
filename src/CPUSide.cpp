/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "System.h"
#include "StaticVals.h"
#include "CPUSide.h" //Spec declaration

CPUSide::CPUSide(System & sys, StaticVals & statV) :
  varRef(sys, statV), pdb(sys, statV), console(varRef), block(varRef),
  hist(varRef), checkpoint(sys, statV)
#if ENSEMBLE == GCMC
  , sample_N_E(varRef)
#endif
{}

void CPUSide::Init(PDBSetup const& pdbSet, config_setup::Output const& out,
                   const ulong tillEquil, const ulong totSteps, ulong startStep)
{
  equilSteps = tillEquil;
  //Initialize arrays in object that collects references and calc'ed vals.
  varRef.Init(pdbSet.atoms);
  //Initialize output components.
  timer.Init(out.console.frequency, totSteps, startStep);
  outObj.push_back(&console);
  outObj.push_back(&pdb);
  if (out.statistics.settings.block.enable)
    outObj.push_back(&block);
  if (out.checkpoint.enable)
    outObj.push_back(&checkpoint);

#if ENSEMBLE == GCMC
  outObj.push_back(&hist);
  outObj.push_back(&sample_N_E);
#endif
  //Calculate pressure, heat of vap. (if applicable), etc.
  varRef.CalcAndConvert(0);
  for (uint o = 0; o < outObj.size(); o++)
    outObj[o]->Init(pdbSet.atoms, out, tillEquil, totSteps);
}

void CPUSide::Output(const ulong step)
{
  //Calculate pressure, heat of vap. (if applicable), etc.
  varRef.CalcAndConvert(step);
  //Do standard output events.
  for (uint o = 0; o < outObj.size(); o++)
    outObj[o]->Output(step);
  timer.CheckTime(step);
}

Clock* CPUSide::getClock(){
  return &timer;
}

void CPUSide::exchangeOfstreamPointers(CPUSide * otherCPUSide){
  ofstream * swapperForConsole = console.getConsoleToFile(); 
  console.setConsoleToFile(otherCPUSide->console.getConsoleToFile());
  otherCPUSide->console.setConsoleToFile(swapperForConsole);

  for (uint box = 0; box < BOX_TOTAL; box++){
    if (pdb.getEnableOutState()){
      ofstream * swapperForPDB = pdb.getPDBToFile(box);
      pdb.setPDBToFile(box, otherCPUSide->pdb.getPDBToFile(box));
      otherCPUSide->pdb.setPDBToFile(box, swapperForPDB);

      if (pdb.getEnableRestOut()){
        ofstream * swapperForPDBRest = pdb.getPDBRestartToFile(box);
        pdb.setPDBRestartToFile(box, otherCPUSide->pdb.getPDBRestartToFile(box));
        otherCPUSide->pdb.setPDBRestartToFile(box, swapperForPDBRest);
      }
    }

    if (block.enableOut){
      ofstream * swapperForBlock = block.getBlockToFile(box);
      block.setBlockToFile(box, otherCPUSide->block.getBlockToFile(box));
      otherCPUSide->block.setBlockToFile(box, swapperForBlock);
    }
#if ENSEMBLE == GCMC

    ofstream * swapperForHist = hist.getHistToFile(box);
    hist.setHistToFile(box, otherCPUSide->hist.getHistToFile(box));
    otherCPUSide->hist.setHistToFile(box, swapperForHist);

    ofstream * swapperForSample_N_E = sample_N_E.getSample_N_EToFile(box);
    sample_N_E.setSample_N_EToFile(box, otherCPUSide->sample_N_E.getSample_N_EToFile(box));
    otherCPUSide->sample_N_E.setSample_N_EToFile(box, swapperForSample_N_E);  

#endif 
  }
}
