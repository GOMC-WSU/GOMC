/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFParticle.h"
#include "NumLib.h" //For Sq, Cb, and MeanA/G functions.
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

FFParticle::FFParticle(Forcefield &ff) : forcefield(ff), mass(NULL), nameFirst(NULL), nameSec(NULL),
  n(NULL), n_1_4(NULL), sigmaSq(NULL), sigmaSq_1_4(NULL), epsilon_cn(NULL),
  epsilon(NULL), epsilon_1_4(NULL), epsilon_cn_1_4(NULL), epsilon_cn_6(NULL),
  epsilon_cn_6_1_4(NULL), nOver6(NULL), nOver6_1_4(NULL), enCorrection(NULL),
  virCorrection(NULL) 
#ifdef GOMC_CUDA
  , varCUDA(NULL)
#endif 
  {
  }

FFParticle::~FFParticle(void)
{
  delete[] mass;
  delete[] nameFirst;
  delete[] nameSec;

  delete[] sigmaSq;
  delete[] n;
  delete[] epsilon;
  delete[] epsilon_cn;
  delete[] epsilon_cn_6;
  delete[] nOver6;
  // parameter for 1-4 interaction, 
  // same one will be used for 1-3 interaction
  delete[] sigmaSq_1_4;
  delete[] n_1_4;
  delete[] epsilon_1_4;
  delete[] epsilon_cn_1_4;
  delete[] epsilon_cn_6_1_4;
  delete[] nOver6_1_4;

  delete[] enCorrection;
  delete[] virCorrection;

#ifdef GOMC_CUDA
  DestroyCUDAVars(varCUDA);
  delete varCUDA;
#endif
}

void FFParticle::Init(ff_setup::Particle const& mie,
                      ff_setup::NBfix const& nbfix)
{
#ifdef GOMC_CUDA
  // Variables for GPU stored in here
  varCUDA = new VariablesCUDA();
#endif
  count = mie.epsilon.size(); //Get # particles read
  //Size LJ particle kind arrays
  mass = new double [count];
  //Size LJ-LJ pair arrays
  uint size = num::Sq(count);
  nameFirst = new std::string [size];
  nameSec = new std::string [size];

#ifdef MIE_INT_ONLY
  n = new uint [size];
  n_1_4 = new uint [size];
#else
  n = new double [size];
  n_1_4 = new double [size];
#endif
  epsilon = new double [size];
  epsilon_cn = new double [size];
  epsilon_cn_6 = new double [size];
  nOver6 = new double [size];
  sigmaSq = new double [size];

  epsilon_1_4 = new double [size];
  epsilon_cn_1_4 = new double [size];
  epsilon_cn_6_1_4 = new double [size];
  nOver6_1_4 = new double [size];
  sigmaSq_1_4 = new double [size];

  enCorrection = new double [size];
  virCorrection = new double [size];

  //Combining VDW parameter
  Blend(mie, forcefield.rCut);
  //Adjusting VDW parameter using NBFIX
  AdjNBfix(mie, nbfix, forcefield.rCut);

#ifdef GOMC_CUDA
  double diElectric_1 = 1.0 / forcefield.dielectric;
  InitGPUForceField(*varCUDA, sigmaSq, epsilon_cn, n, forcefield.vdwKind,
                    forcefield.isMartini, count, forcefield.rCut,
		                forcefield.rCutCoulomb,forcefield.rCutLow,
		                forcefield.rswitch, forcefield.alpha, forcefield.ewald,
		                diElectric_1);
#endif
}

double FFParticle::EnergyLRC(const uint kind1, const uint kind2) const
{
  return enCorrection[FlatIndex(kind1, kind2)];
}

double FFParticle::VirialLRC(const uint kind1, const uint kind2) const
{
  return virCorrection[FlatIndex(kind1, kind2)];
}

void FFParticle::Blend(ff_setup::Particle const& mie, const double rCut)
{
  for(uint i = 0; i < count; ++i) {
    for(uint j = 0; j < count; ++j) {
      uint idx = FlatIndex(i, j);
      //get all name combination for using in nbfix
      nameFirst[idx] = mie.getname(i);
      nameFirst[idx] += mie.getname(j);
      nameSec[idx] = mie.getname(j);
      nameSec[idx] += mie.getname(i);

      n[idx] = num::MeanA(mie.n, mie.n, i, j);
      n_1_4[idx] = num::MeanA(mie.n_1_4, mie.n_1_4, i, j);
      double cn = n[idx] / (n[idx] - 6) * pow(n[idx] / 6, (6 / (n[idx] - 6)));
      double cn_1_4 = n_1_4[idx] / (n_1_4[idx] - 6) *
                      pow(n_1_4[idx] / 6, (6 / (n_1_4[idx] - 6)));

      double sigma, sigma_1_4;
      sigma = sigma_1_4 = 0.0;
      
      if(forcefield.vdwGeometricSigma) {
        sigma = num::MeanG(mie.sigma, mie.sigma, i, j);
        sigma_1_4 = num::MeanG(mie.sigma_1_4, mie.sigma_1_4, i, j);
      } else {
        sigma = num::MeanA(mie.sigma, mie.sigma, i, j);
        sigma_1_4 = num::MeanA(mie.sigma_1_4, mie.sigma_1_4, i, j);
      }

      double tc = 1.0;
      double rRat = sigma / rCut;
      // calculate sig^2 and tc*sig^3
      num::Cb(sigmaSq[idx], tc, sigma);
      sigmaSq_1_4[idx] = sigma_1_4 * sigma_1_4;
      tc *= 0.5 * 4.0 * M_PI;
      epsilon[idx] = num::MeanG(mie.epsilon, mie.epsilon, i, j);
      epsilon_cn[idx] = cn * epsilon[idx];
      epsilon_1_4[idx] = num::MeanG(mie.epsilon_1_4, mie.epsilon_1_4, i, j);
      epsilon_cn_1_4[idx] = cn * epsilon_1_4[idx];
      epsilon_cn_6[idx] = epsilon_cn[idx] * 6;
      epsilon_cn_6_1_4[idx] = epsilon_cn_1_4[idx] * 6;
      nOver6[idx] = n[idx] / 6;
      nOver6_1_4[idx] = n_1_4[idx] / 6;
      enCorrection[idx] = tc / (n[idx] - 3) * epsilon_cn[idx] *
                          ( pow(rRat, n[idx] - 3) -
                            (double)(n[idx] - 3.0) / 3.0 * pow(rRat, 3) );
      virCorrection[idx] = tc / (n[idx] - 3) * epsilon_cn_6[idx] *
                           ( (double)(n[idx]) / 6.0 * pow(rRat, n[idx] - 3) -
                             (double)(n[idx] - 3.0) / 3.0 * pow(rRat, 3) );
    }
  }
}

void FFParticle::AdjNBfix(ff_setup::Particle const& mie,
                          ff_setup::NBfix const& nbfix, const double rCut)
{
  uint size = num::Sq(count);
  for(uint i = 0; i < nbfix.epsilon.size(); i++) {
    for(uint j = 0; j < size; j++) {
      if(nbfix.getname(i) == nameFirst[j] || nbfix.getname(i) ==  nameSec[j]) {
        n[j] = nbfix.n[i];
        n_1_4[j] = nbfix.n_1_4[i];
        double rRat = nbfix.sigma[i] / rCut, tc = 1.0;
        //calculating sig^2 and tc*sig^3
        num::Cb(sigmaSq[j], tc, nbfix.sigma[i]);
        sigmaSq_1_4[j] = nbfix.sigma_1_4[i] * nbfix.sigma_1_4[i];
        tc *= 0.5 * 4.0 * M_PI;
        double cn = n[j] / (n[j] - 6) * pow(n[j] / 6, (6 / (n[j] - 6)));
        double cn_1_4 = n_1_4[j] / (n_1_4[j] - 6) *
                        pow(n_1_4[j] / 6, (6 / (n_1_4[j] - 6)));
        epsilon[j] = nbfix.epsilon[i];
        epsilon_cn[j] = cn * nbfix.epsilon[i];
        epsilon_1_4[j] = nbfix.epsilon_1_4[i];
        epsilon_cn_1_4[j] = cn_1_4 * nbfix.epsilon_1_4[i];
        epsilon_cn_6[j] = epsilon_cn[j] * 6;
        epsilon_cn_6_1_4[j] = epsilon_cn_1_4[j] * 6;
        nOver6[j] = n[j] / 6;
        nOver6_1_4[j] = n_1_4[j] / 6;
        enCorrection[j] = tc / (n[j] - 3) * epsilon_cn[j] *
                          ( pow(rRat, n[j] - 3) -
                            (double)(n[j] - 3.0) / 3.0 * pow(rRat, 3) );
        virCorrection[j] = tc / (n[j] - 3) * epsilon_cn_6[j] *
                           ( (double)(n[j]) / 6.0 * pow(rRat, n[j] - 3) -
                             (double)(n[j] - 3.0) / 3.0 * pow(rRat, 3) );
      }
    }
  }
}

double FFParticle::GetEpsilon(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return epsilon[idx];
}
double FFParticle::GetEpsilon_1_4(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return epsilon_1_4[idx];
}
double FFParticle::GetSigma(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return sqrt(sigmaSq[idx]);
}
double FFParticle::GetSigma_1_4(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return sqrt(sigmaSq_1_4[idx]);
}
double FFParticle::GetN(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return n[idx];
}
double FFParticle::GetN_1_4(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return n_1_4[idx];
}

// Defining the functions

inline void FFParticle::CalcAdd_1_4(double& en, const double distSq,
                                    const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq_1_4[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n_1_4[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n_1_4[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  en += epsilon_cn_1_4[index] * (repulse - attract);
}

inline void FFParticle::CalcCoulombAdd_1_4(double& en, const double distSq,
    const double qi_qj_Fact,
    const bool NB) const
{
  double dist = sqrt(distSq);
  if(NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}



//mie potential
inline double FFParticle::CalcEn(const double distSq,
                                 const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  return epsilon_cn[index] * (repulse - attract);
}

double FFParticle::CalcEnAttract(const double distSq) const
{
  if (forcefield.rCutSq < distSq)
    return 0.0;
  return 1.0 / (distSq * distSq * distSq);
}

double FFParticle::CalcEnRepulse(const double distSq) const
{
  if (forcefield.rCutSq < distSq)
    return 0.0;
  double r6 = 1.0 / (distSq * distSq * distSq);
  return r6 * r6;
}

inline double FFParticle::CalcCoulomb(const double distSq,
                                      const double qi_qj_Fact, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return  qi_qj_Fact * erfc(val) / dist;
  } else {
    double dist = sqrt(distSq);
    return  qi_qj_Fact / dist;
  }
}

double FFParticle::CalcCoulombNoFact(const double distSq, const uint b) const
{
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if (forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return  erfc(val) / dist;
  }
  else {
    return sqrt(distSq);
  }
}

inline double FFParticle::CalcVir(const double distSq,
                                  const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double rNeg2 = 1.0 / distSq;
  double rRat2 = rNeg2 * sigmaSq[index];
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  //Virial is the derivative of the pressure... mu
  return epsilon_cn_6[index] * (nOver6[index] * repulse - attract) * rNeg2;
}

inline double FFParticle::CalcCoulombVir(const double distSq,
    const double qi_qj, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double constValue = 2.0 * forcefield.alpha[b] / sqrt(M_PI);
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = 1.0 - erf(forcefield.alpha[b] * dist);
    return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    double dist = sqrt(distSq);
    return qi_qj / (distSq * dist);
  }
}

inline double FFParticle::CalcCoulombVirNoFact(const double distSq,
  const uint b) const
{
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if (forcefield.ewald) {
    double dist = sqrt(distSq);
    double constValue = 2.0 * forcefield.alpha[b] / sqrt(M_PI);
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = 1.0 - erf(forcefield.alpha[b] * dist);
    return (temp / dist + constValue * expConstValue) / distSq;
  }
  else {
    double dist = sqrt(distSq);
    return 1.0 / (distSq * dist);
  }
}

void FFParticle::InitializeTables()
{
  double r;
  energyTable.resize(BOXES_WITH_U_NB);

  for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
    // Initialize loca variables
    int tableLength           = ((forcefield.rCutSq + 1) / TABLE_STEP);
    uint kindTotal            = count;
    uint kindTotalSq          = kindTotal * kindTotal;

    // Allocate space for temporary ar23rays
    LJAttractV.resize(tableLength);
    LJAttractF.resize(tableLength);
    LJRepulseV.resize(tableLength);
    LJRepulseF.resize(tableLength);
    CoulombV.resize(tableLength);
    CoulombF.resize(tableLength);

    // Let's make sure we are not allocating too much data for energy table.
    if (tableLength * kindTotalSq > 10000000) {
      std::cerr << "Error: There is more than 1 million entries in tabulated potential.\n";
      std::cerr << "       Please turn off energy tables.\n";
      exit(EXIT_FAILURE);
    }

    // Allocate tables for energy and force    
    energyTable[b].resize(kindTotalSq * tableLength * TABLE_STRIDE);
    
    for (uint k1 = 0; k1 < kindTotal; k1++) {
      for (uint k2 = 0; k2 < kindTotal; k2++) {
        uint k = FlatIndex(k1, k2);
        for (int i = 0; i < tableLength; i++) {
          r = i * TABLE_STEP;

          // Firt 4 indeces are for Attraction of LJ
          LJAttractV[i] = CalcEnAttract(r * r);
          LJAttractF[i] = LJAttractV[i] * 6.0 / r;
          LJRepulseV[i] = CalcEnRepulse(r * r);
          LJRepulseF[i] = LJRepulseV[i] * n[k] / r;
          CoulombV[i] = CalcCoulombNoFact(r * r, b);
          CoulombF[i] = -1 * CalcCoulombNoFact(r * r, b) / r;
        }

        for (int i = 0; i < tableLength; i++) {
          int index = k * tableLength + i * TABLE_STRIDE;
          r = i * TABLE_STEP;

          if (i < tableLength - 1) {
            energyTable[b][index]      = TABLE_STEP;
            energyTable[b][index + 1]  = -1.0 * LJAttractF[i] * TABLE_STEP;
            energyTable[b][index + 2]  = 3 * (LJAttractV[i + 1] - LJAttractV[i]) +
                                         (LJAttractF[i + 1] - 2 * LJAttractF[i]) *
                                         TABLE_STEP;
            energyTable[b][index + 3]  = -2 * (LJAttractV[i + 1] - LJAttractV[i]) -
                                         (LJAttractF[i + 1] + LJAttractF[i]) *
                                         TABLE_STEP;
            energyTable[b][index + 4]  = TABLE_STEP;
            energyTable[b][index + 5]  = -1.0 * LJRepulseF[i] * TABLE_STEP;
            energyTable[b][index + 6]  = 3 * (LJRepulseV[i + 1] - LJRepulseV[i]) +
                                         (LJRepulseF[i + 1] - 2 * LJRepulseF[i]) *
                                         TABLE_STEP;
            energyTable[b][index + 7]  = -2 * (LJRepulseV[i + 1] - LJRepulseV[i]) -
                                         (LJRepulseF[i + 1] + LJRepulseF[i]) *
                                         TABLE_STEP;
            energyTable[b][index + 8]  = TABLE_STEP;
            energyTable[b][index + 9]  = -1.0 * CoulombF[i] * TABLE_STEP;
            energyTable[b][index + 10] = 3 * (CoulombV[i + 1] - CoulombV[i]) +
                                         (CoulombF[i + 1] - 2 * CoulombF[i]) *
                                         TABLE_STEP;
            energyTable[b][index + 11] = -2 * (CoulombV[i + 1] - CoulombV[i]) -
                                         (CoulombF[i + 1] + CoulombF[i]) *
                                         TABLE_STEP;
          }
          else {
            energyTable[b][index]      = TABLE_STEP;
            energyTable[b][index + 1]  = -1.0 * LJAttractF[i] * TABLE_STEP;
            energyTable[b][index + 2]  = 0.0;
            energyTable[b][index + 3]  = 0.0;
            energyTable[b][index + 4]  = TABLE_STEP;
            energyTable[b][index + 5]  = -1.0 * LJRepulseF[i] * TABLE_STEP;
            energyTable[b][index + 6]  = 0.0;
            energyTable[b][index + 7]  = 0.0;
            energyTable[b][index + 8]  = TABLE_STEP;
            energyTable[b][index + 9]  = -1.0 * CoulombF[i] * TABLE_STEP;
            energyTable[b][index + 10] = 0.0;
            energyTable[b][index + 11] = 0.0;
          }
        }
      }
    }
  }
}

double FFParticle::CalcEnEnergyTable(const double distSq,
                                     const uint kind1, const uint kind2) const
{

}