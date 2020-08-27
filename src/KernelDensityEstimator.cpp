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
  for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
    outF[b] = NULL;
    name[b] = NULL;
  }

  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    energiesPerBox[b].resize(sampler.samplesPerFrame);
    probabilitiesPerBox[b].resize(sampler.samplesPerFrame);
  }
}

KernelDensityEstimator::~KernelDensityEstimator()
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b]->is_open()) outF[b]->close();
  }
}

void KernelDensityEstimator::GeneratePDF(){

  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    std::copy(sampler.samplesE[b], sampler.samplesE[b] + sampler.samplesPerFrame, energiesPerBox[b].begin());
    std::sort(energiesPerBox[b].begin(), energiesPerBox[b].end());
  }

  switch(k_h){
    case BOXCAR:
      int i, b, myBinCount, distance;
      #if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
      #pragma omp parallel for default(none) private(b, i, myBinCount, distance) shared(energiesPerBox, probabilitiesPerBox) collapse(2)
      #endif
      for(b = 0; b < BOXES_WITH_U_NB; b++){
        for (i = 0; i < sampler.samplesPerFrame; i++){
          myBinCount = 1;
          distance = 1;
          while (i < sampler.samplesPerFrame && i + distance < sampler.samplesPerFrame &&  abs(energiesPerBox[b][i] - energiesPerBox[b][i + distance]) < h ){
            myBinCount++;
            distance++;
          }
          distance = 1;
          while (i > 0 && i - distance >= 0 &&  abs(energiesPerBox[b][i] - energiesPerBox[b][i - distance]) < h ){
            myBinCount++;
            distance++;
          }
          probabilitiesPerBox[b][i] =  (double)myBinCount / (double)sampler.samplesPerFrame;
        }
      }

      break;

    case GAUSSIAN:

      break;

    default:

      break;
  }
}
