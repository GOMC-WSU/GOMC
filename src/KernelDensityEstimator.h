/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H
  
#include <string>
#include <fstream>
#include "EnPartCntSampleOutput.h"

enum KernelType{ BOXCAR, GAUSSIAN };


class KernelDensityEstimator : public OutputableBase {
  public:
  KernelDensityEstimator(EnPartCntSample & sampler);
  void GeneratePDF();
  ~KernelDensityEstimator();

  //We use the EnPartCntSample as our sampler, so do nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

  virtual void DoOutput(const ulong step);

  void PrintKDE(const uint b);  

  
  uint kdeCounter;
  double h;
  KernelType k_h;
  uint stepsPerSample;
  std::ofstream outF [BOXES_WITH_U_NB];
  std::string name [BOXES_WITH_U_NB];
  EnPartCntSample & sampler;
  std::vector< std::vector<double> > energiesPerBox;
  std::vector< std::vector<double> > probabilitiesPerBox;
  friend class EnPartCntSample;
};

#endif /*KERNELDENSITYESTIMATOR_H*/