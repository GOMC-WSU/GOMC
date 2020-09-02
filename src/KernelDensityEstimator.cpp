/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "KernelDensityEstimator.h"
#include <sstream>

KernelDensityEstimator::KernelDensityEstimator(EnPartCntSample & enPCS) : sampler(enPCS)
{
  kdeCounter = 0;
  h = 1.0;
}

KernelDensityEstimator::~KernelDensityEstimator()
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open()) outF[b].close();
  }
}

void KernelDensityEstimator::Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output)
{
  energiesPerBox.resize(BOXES_WITH_U_NB);
  probabilitiesPerBox.resize(BOXES_WITH_U_NB);
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    energiesPerBox[b].resize(sampler.samplesPerFrame);
    probabilitiesPerBox[b].resize(sampler.samplesPerFrame);
  }
  std::string bStr = "", aliasStr = "", numStr = "";
  sstrm::Converter toStr;
  stepsPerSample = output.state.files.hist.stepsPerHistSample;
  stepsPerOut = output.statistics.settings.hist.frequency;
  enableOut = output.statistics.settings.hist.enable;
  if (enableOut) {
    //Assign arrays for boxes of interest
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      //Get alias string, based on box #.
      bStr = "Box ";
      numStr = "";
      toStr << b + 1;
      toStr >> numStr;
      aliasStr = "Output KDE file for Box ";
      aliasStr += numStr;
      bool notify;
#ifndef NDEBUG
      notify = true;
#else
      notify = false;
#endif
#if GOMC_LIB_MPI
      name[b] = pathToReplicaDirectory + "kde_box_" + std::to_string(b);
#else
      name[b] = "kde_box_" + std::to_string(b);
#endif
      outF[b].open(name[b].c_str(), std::ofstream::out);
    } 
  }
}

void KernelDensityEstimator::GeneratePDF(){

  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    std::copy(sampler.samplesE[b], sampler.samplesE[b] + sampler.kdeSamplesCollectedInFrame, energiesPerBox[b].begin());
    std::sort(energiesPerBox[b].begin(), energiesPerBox[b].end());
  }

  switch(k_h){
    case BOXCAR:
      int i, b, myBinCount, distance;
      #if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
      #pragma omp parallel for default(none) private(b, i, myBinCount, distance) shared(energiesPerBox, probabilitiesPerBox) collapse(2)
      #endif
      for(b = 0; b < BOXES_WITH_U_NB; b++){
        for (i = 0; i < sampler.kdeSamplesCollectedInFrame; i++){
          myBinCount = 1;
          distance = 1;
          while (i < sampler.kdeSamplesCollectedInFrame && i + distance < sampler.kdeSamplesCollectedInFrame &&  abs(energiesPerBox[b][i] - energiesPerBox[b][i + distance]) < h ){
            myBinCount++;
            distance++;
          }
          distance = 1;
          while (i > 0 && i - distance >= 0 &&  abs(energiesPerBox[b][i] - energiesPerBox[b][i - distance]) < h ){
            myBinCount++;
            distance++;
          }
          probabilitiesPerBox[b][i] =  (double)myBinCount / (double)sampler.kdeSamplesCollectedInFrame;
        }
      }

      break;

    case GAUSSIAN:

      break;

    default:

      break;
  }
}

void KernelDensityEstimator::DoOutput(const ulong step)
{
  //Don't output until equilibrated.
  if ((step) < stepsTillEquil) return;
  //Write to histogram file, if equilibrated.
  if ((step + 1) % stepsPerOut == 0) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
        if (outF[b].is_open()){
          GeneratePDF();
          PrintKDE(b);
        } else{
          std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                    << "(kde file)" << std::endl;
          outF[b].close();
        }
    }
  }
}


void KernelDensityEstimator::PrintKDE(const uint b)
{
    for (uint n = 0; n < sampler.kdeSamplesCollectedInFrame; ++n) {
        outF[b].precision(dbl::max_digits10);
        outF[b] << kdeCounter << " " << energiesPerBox[b][n] << " " << probabilitiesPerBox[b][n] << std::endl;
    }
    kdeCounter++;
}