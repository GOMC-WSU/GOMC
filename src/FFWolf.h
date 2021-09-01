/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef WOLF_H
#define WOLF_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "Forcefield.h"
#include "BasicTypes.h" //for uint
#include "NumLib.h" //For Cb, Sq
#include "Setup.h"
#ifdef GOMC_CUDA
#include "VariablesCUDA.cuh"
#endif


//////////////////////////////////////////////////////////////////////
////////////////////////// Wolf Style ////////////////////////////
//////////////////////////////////////////////////////////////////////
// LJ potential calculation:
// U_LJ(rij) = 4 * eps_ij * ( (sig_ij/rij)^12 - (sig_ij/rij)^6)
// U_LJ_tc = 1/2 * sum i,j [(16 * pi * N_i * N_j * eps_ij / V) * (((sig_ij^12)/9*R_c^9) - ((sig_ij^6)/3*R_c^3))]
//
// U_elect = (1) - (2)
// (1) 1 / (4*pi*eps_0) * sum_i<j,_rij<Rc [qi * qj ((erfc(alpha*r_ij)/r_ij) - (erfc(alpha*Rc)/Rc))] (1)
// (2) 1 / (4*pi*eps_0) * (erfc(alpha*Rc)/2Rc + alpha/sqrt(pi)) * sum i = 1,N [(q_i) ^2]


struct FF_WOLF {
public:
  FF_WOLF(Forcefield &ff)
  {
    alpha = lambda = thermalWavelength = 0.0;
  }
protected:
    double alpha, lambda, thermalWavelength;
    virtual double CalcEn(const double distSq,
                          const uint kind1, const uint kind2,
                          const double lambda) const;
    virtual double CalcTailCorrection(const double distSq, const uint index) const;
    virtual double CalcCoulomb(const double distSq, const double qi_qj_Fact,
                              const uint b) const;


    double *shiftConst, *shiftConst_1_4;
};

#endif