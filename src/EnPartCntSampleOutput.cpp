/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "EnPartCntSampleOutput.h"

#include <limits> // for std::numeric_limit

#include "ConfigSetup.h"
#include "OutConst.h"
#include "OutputVars.h"
#include "PDBConst.h"
#include "System.h"

#if ENSEMBLE == GCMC

EnPartCntSample::~EnPartCntSample() {
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open())
      outF[b].close();
    if (samplesE[b] != NULL)
      delete[] samplesE[b];
    if (samplesN[b] != NULL) {
      for (uint k = 0; k < var->numKinds; ++k)
        if (samplesN[b][k] != NULL)
          delete[] samplesN[b][k];
    }
    if (samplesN[b] != NULL)
      delete[] samplesN[b];
  }
}

void EnPartCntSample::Init(pdb_setup::Atoms const &atoms,
                           config_setup::Output const &output) {
  InitVals(output.statistics.settings.hist);
  if (enableOut) {
    stepsPerSample = output.state.files.hist.stepsPerHistSample;
    uint samplesPerFrame =
        output.statistics.settings.hist.frequency / stepsPerSample + 1;
    samplesCollectedInFrame = 0;
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
#if GOMC_LIB_MPI
      name[b] = pathToReplicaOutputDirectory +
                GetFName(output.state.files.hist.sampleName,
                         output.state.files.hist.number,
                         output.state.files.hist.letter, b);
#else
      name[b] = GetFName(output.state.files.hist.sampleName,
                         output.state.files.hist.number,
                         output.state.files.hist.letter, b);
#endif
      samplesE[b] = new double[samplesPerFrame];
      samplesN[b] = new uint *[var->numKinds];
      for (uint k = 0; k < var->numKinds; ++k) {
        samplesN[b][k] = new uint[samplesPerFrame];
      }
      outF[b].open(name[b].c_str(), std::ofstream::out);
    }
    WriteHeader();
  }
}

void EnPartCntSample::Sample(const ulong step) {
  // Don't sample until equilibrated.
  if ((step) < stepsTillEquil)
    return;
  // Only sample on specified interval.
  if ((step + 1) % stepsPerSample == 0) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      samplesE[b][samplesCollectedInFrame] = var->energyRef[b].inter +
                                             var->energyRef[b].tailCorrection +
                                             var->energyRef[b].totalElect;

      for (uint k = 0; k < var->numKinds; ++k) {
        samplesN[b][k][samplesCollectedInFrame] =
            var->numByKindBox[k + var->numKinds * b];
      }
    }
    ++samplesCollectedInFrame;
  }
}

void EnPartCntSample::WriteHeader(void) {
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open()) {
      outF[b] << var->T_in_K << " " << var->numKinds << " ";
#if ENSEMBLE == GCMC
      for (uint k = 0; k < var->numKinds; k++) {
        outF[b] << var->kindsRef[k].chemPot << " ";
      }
#endif
      XYZ bAx = var->axisRef->Get(0);
      outF[b] << bAx.x << " " << bAx.y << " " << bAx.z << std::endl;
      outF[b] << std::setprecision(std::numeric_limits<double>::digits10 + 2);
      outF[b].setf(std::ios_base::left, std::ios_base::adjustfield);
    } else
      std::cerr << "Unable to write to file \"" << name[b] << "\" "
                << "(energy and part. num samples file)" << std::endl;
  }
}

void EnPartCntSample::DoOutput(const ulong step) {
  // Don't output until equilibrated.  Also don't output if this we continue
  // a checkpoint or init directly on an output step.
  if ((step) < stepsTillEquil ||
      ((restartFromCheckpoint || initStepRead) && step == startStep))
    return;
  // Output a sample in the form <N1,... Nk, E_total>
  // Only sample on specified interval.
  if ((step + 1) % stepsPerOut == 0) {
    GOMC_EVENT_START(1, GomcProfileEvent::HIST_OUTPUT);
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      if (outF[b].is_open()) {
        for (uint n = 0; n < samplesCollectedInFrame; ++n) {
          for (uint k = 0; k < var->numKinds; k++) {
            outF[b] << std::setw(11) << samplesN[b][k][n] << " ";
          }
          outF[b] << std::setw(25) << samplesE[b][n] << std::endl;
        }
      } else
        std::cerr << "Unable to write to file \"" << name[b] << "\" "
                  << "(energy and part. num samples file)" << std::endl;
    }
    GOMC_EVENT_STOP(1, GomcProfileEvent::HIST_OUTPUT);
  }
  samplesCollectedInFrame = 0;
}

void EnPartCntSample::DoOutputRestart(const ulong step) {}

std::string EnPartCntSample::GetFName(std::string const &sampleName,
                                      std::string const &histNum,
                                      std::string const &histLetter,
                                      const uint box) {
  std::stringstream sstrm;
  std::string fName = sampleName, strBox;
  fName += histNum;
  fName += histLetter;
  if (BOXES_WITH_U_NB > 1) {
    fName += "_box";
    sstrm << box;
    sstrm >> strBox;
    fName += strBox;
  }
  fName += ".dat";
  return fName;
}

#endif /*ENSEMBLE==GCMC*/
