/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "EnPartCntSampleOutput.h"
#include "PDBConst.h"
#include "OutConst.h"
#include "OutputVars.h"
#include "ConfigSetup.h"
#include "System.h"

#if ENSEMBLE == GCMC

EnPartCntSample::~EnPartCntSample()
{
   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      if (samplesE[b] != NULL) delete[] samplesE[b];
      for (uint k = 0; k < var->numKinds; ++k)
         if (samplesN[b][k] != NULL) delete[] samplesN[b][k];
      if (samplesN[b] != NULL) delete[] samplesN[b];
   }
}

void EnPartCntSample::Init(pdb_setup::Atoms const& atoms,
                           config_setup::Output const& output)
{
   InitVals(output.statistics.settings.hist);
   if (enableOut)
   {
      stepsPerSample = output.statistics.settings.hist.frequency /
         output.state.files.hist.samplesPerHist;
      samplesCollectedInFrame = 0;
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
      {
         name[b] = GetFName(output.state.files.hist.sampleName,
                            output.state.files.hist.number,
                            output.state.files.hist.letter,
                            b);
         samplesE[b] = new double [output.state.files.hist.samplesPerHist+1];
         samplesN[b] = new uint * [var->numKinds];
         for (uint k = 0; k < var->numKinds; ++k)
         {
            samplesN[b][k] =
               new uint [output.state.files.hist.samplesPerHist+1];
         }
      }
      WriteHeader();
   }
}

void EnPartCntSample::Sample(const ulong step)
{
   //Don't sample until equilibrated.
   if ((step+1) < stepsTillEquil) return;
   //Only sample on specified interval.
   if ((step+1) % stepsPerSample == 0)
   {
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
      {
         samplesE[b][samplesCollectedInFrame] = var->energyRef[b].total;
         for (uint k = 0; k < var->numKinds; ++k)
         {
            samplesN[b][k][samplesCollectedInFrame] =
               var->numByKindBox[k+var->numKinds*b];
         }
      }
      ++samplesCollectedInFrame;
   }
}

void EnPartCntSample::WriteHeader(void)
{
   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      outF.open(name[b].c_str(), std::ofstream::out);
      if (outF.is_open())
      {
         outF << var->T_in_K << " " << var->numKinds << " ";
#if ENSEMBLE == GCMC
         for (uint k = 0; k < var->numKinds; k++)
         {
            outF << var->kindsRef[k].chemPot << " ";
         }
#endif
         XYZ bAx = var->axisRef->Get(0);
         outF << bAx.x << " " << bAx.y << " " << bAx.z << std::endl;
         outF.close();
      }
      else
         std::cerr << "Unable to write to file \"" <<  name[b] << "\" " 
                   << "(energy and part. num samples file)" << std::endl;
   }
}

void EnPartCntSample::DoOutput(const ulong step)
{
   //Don't output until equilibrated.
   if ((step+1) < stepsTillEquil) return;
   //Output a sample in the form <N1,... Nk, E_total>
   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      outF.open(name[b].c_str(), std::ofstream::app);
      if (outF.is_open())
      {
         for (uint n = 0; n < samplesCollectedInFrame; ++n)
         {
            for (uint k = 0; k < var->numKinds; k++)
            {
               outF << samplesN[b][k][n] << " ";
            }
            outF << samplesE[b][n] << std::endl;
         }
      }
      else
         std::cerr << "Unable to write to file \"" <<  name[b] << "\" " 
                   << "(energy and part. num samples file)" << std::endl;
      outF.close();
   }
   samplesCollectedInFrame = 0;
}

   
std::string EnPartCntSample::GetFName(std::string const& sampleName, 
                                      std::string const& histNum,
                                      std::string const& histLetter,
                                      const uint box)
{ 
   std::stringstream sstrm;
   std::string fName = sampleName, strBox;
   fName += histNum;
   fName += histLetter;
   if ( BOXES_WITH_U_NB > 1 )
   {
      fName += "_box";
      sstrm << box;
      sstrm >> strBox;
      fName += strBox;
   }
   fName += ".dat";
   return fName;
}

#endif /*ENSEMBLE==GCMC*/

