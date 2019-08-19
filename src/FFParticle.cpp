/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFParticle.h"
#include "NumLib.h" //For Sq, Cb, and MeanA/G functions.
#include <functional> // for bind2nd
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

FFParticle::FFParticle(Forcefield &ff) : forcefield(ff), mass(NULL), nameFirst(NULL), nameSec(NULL),
  n(NULL), n_1_4(NULL), sigmaSq(NULL), epsilon(NULL), epsilon_cn(NULL), epsilon_1_4(NULL),
  sigmaSq_1_4(NULL), epsilon_cn_1_4(NULL), epsilon_cn_6(NULL),
  epsilon_cn_6_1_4(NULL), nOver6(NULL), nOver6_1_4(NULL), enCorrection(NULL),
  virCorrection(NULL) 
#ifdef GOMC_CUDA
  , varCUDA(NULL)
#endif 
{
  energyTableEnabled = false;
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
  if(energyTableEnabled) {
    return CSTable_CalcEnAttract[mv::BOX0][index](distSq) +
      CSTable_CalcEnRepulse[mv::BOX0][index](distSq);
  }

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

inline double FFParticle::CalcEnAttract(const double distSq,
                                           const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;

  return epsilon_cn[index] * -attract;
}

inline double FFParticle::CalcEnRepulse(const double distSq,
                                 const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq[index] / distSq;
#ifdef MIE_INT_ONLY
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  uint n_ij = n[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  return epsilon_cn[index] * repulse;
}

inline double FFParticle::CalcCoulomb(const double distSq,
                                      const double qi_qj_Fact, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;
  
  if(energyTableEnabled) {
    return qi_qj_Fact * CSTable_CalcCoulomb[b][0](distSq);
  }

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
    double dist = sqrt(distSq);
    return 1.0 / dist;
  }
}

// derivative of CalcCoulomb without qqfact 
// ewald: -erfc(a*x)/x^2-(2*a*e^(-a^2*x^2))/(sqrt(pi)*x)
// noewa: -1/x^2
double FFParticle::CalcCoulombNoFactDerivative(const double distSq, const uint b) const
{
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  double a = forcefield.alpha[b];
  if (forcefield.ewald) {
    double x = sqrt(distSq);
    double first = erfc(a * x) / distSq;
    double second_top = 2 * a * exp(-1.0 * pow(a, 2.0) * pow(x, 2.0));
    double second_bot = sqrt(M_PI) * x;
    double second = second_top / second_bot;
    return -(first+second);
  }
  else {
    return -1.0 / distSq;
  }
}

inline double FFParticle::CalcVir(const double distSq,
                                  const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);

  if(energyTableEnabled) {
    return CSTable_CalcEnAttract[mv::BOX0][index].GetDerivativeValue(distSq) +
      CSTable_CalcEnRepulse[mv::BOX0][index].GetDerivativeValue(distSq);
  }

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

  if(energyTableEnabled) {
    return qi_qj * CSTable_CalcCoulomb[b][0].GetDerivativeValue(distSq);
  }

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
    double temp = erfc(forcefield.alpha[b] * dist);
    return (temp / dist + constValue * expConstValue) / distSq;
    // y = (erfc(a*x)/x+e^(-a*x^2)+(2*a)/sqrt(pi))/x^2
  }
  else {
    double dist = sqrt(distSq);
    return 1.0 / (distSq * dist);
    // y = 1.0 / x^3
  }
}

// derivative of above expression. A little complicated!
// ewald: (-((2 x (2 a+a E^(-a^2 x^2)+E^(-a x^2) Sqrt[\[Pi]] (1+a x^2)))/Sqrt[\[Pi]])-3 Erfc[a x])/x^4
// noewa: -3/x^4
double FFParticle::CalcCoulombVirNoFactDerivative(const double distSq,
  const uint b) const
{
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  double dist = sqrt(distSq);
  double a = forcefield.alpha[b];
  double a2 = a * a;
  double sp = sqrt(M_PI);
  if (forcefield.ewald) {
    double temp1 = a * exp(-1.0 * a2 * distSq);
    double temp2 = sp * exp(-1.0 * a * distSq);
    double temp3 = a * distSq + 1;
    double temp4 = 2 * dist * (temp1 + temp2 * temp3 + 2.0 * a);
    double top   = -(temp4 / sp) - 3.0 * erfc(a * dist);
    double ans   = top / (distSq * distSq);
    return ans;
  }
  else {
    return -3.0 / (distSq * distSq);
  }
}

void FFParticle::InitializeTables()
{
  energyTableEnabled = true;

  double r2, r1;
  CSTable_CalcCoulomb      = new CubicSpline*[BOXES_WITH_U_B];
  CSTable_CalcEnAttract    = new CubicSpline*[BOXES_WITH_U_B];
  CSTable_CalcEnRepulse    = new CubicSpline*[BOXES_WITH_U_B];

  for (uint b = 0; b < BOXES_WITH_U_B; b++) {
    // Initialize loca variables
    int tableLength           = ((forcefield.rCutSq + 1) / TABLE_STEP);
    uint kindTotal            = count;
    uint kindTotalSq          = kindTotal * kindTotal;
    CSTable_CalcCoulomb[b]    = new CubicSpline[kindTotalSq];
    CSTable_CalcEnAttract[b]  = new CubicSpline[kindTotalSq];
    CSTable_CalcEnRepulse[b]  = new CubicSpline[kindTotalSq];

    // Let's make sure we are not allocating too much data for energy table.
    if (tableLength * kindTotalSq > 10000000) {
      std::cerr << "Error: There is more than 1 million entries in tabulated potential.\n";
      std::cerr << "       Please turn off energy tables.\n";
      exit(EXIT_FAILURE);
    }
    
    for (uint k1 = 0; k1 < kindTotal; k1++) {
      for (uint k2 = 0; k2 < kindTotal; k2++) {
        // initialize variables
        uint k = FlatIndex(k1, k2);

        CSTable_CalcCoulomb[b][k].Reconstruct(tableLength, TABLE_STEP, 0.0);
        CSTable_CalcEnAttract[b][k].Reconstruct(tableLength, TABLE_STEP, 0.0);
        CSTable_CalcEnRepulse[b][k].Reconstruct(tableLength, TABLE_STEP, 0.0);

        std::vector<double> v_coulomb, v_d_coulomb, v_attract, v_d_attract, v_repulse, v_d_repulse;
        v_coulomb.resize(tableLength);
        v_d_coulomb.resize(tableLength);
        v_attract.resize(tableLength);
        v_d_attract.resize(tableLength);
        v_repulse.resize(tableLength);
        v_d_repulse.resize(tableLength);

        for (int i = 0; i < tableLength; i++) {
          r2 = i * TABLE_STEP;
          r1 = sqrt(r2);

          if(r1 < 0.4) {
            v_coulomb[i]         = 0.0;
            v_d_coulomb[i]       = 0.0;
            v_attract[i]         = 0.0;
            v_d_attract[i]       = 0.0;
            v_repulse[i]         = 0.0;
            v_d_repulse[i]       = 0.0;
          }

          // Firt 4 indeces are for Attraction of LJ
          v_coulomb[i]         = CalcCoulombNoFact(r2, b);
          v_d_coulomb[i]       = CalcCoulombNoFactDerivative(r2, b);
          v_attract[i]         = CalcEnAttract(r2, k1, k2);
          v_d_attract[i]       = -6.0 * v_attract[i] / r1;
          v_repulse[i]         = CalcEnRepulse(r2, k1, k2);
          v_d_repulse[i]       = -n[k] * v_repulse[i] / r1;
        }

        for (int i = 0; i < tableLength; i++) {
          if (i < tableLength - 1) {
            CSTable_CalcCoulomb[b][k].InitializeSpecificPoint(v_coulomb[i], v_coulomb[i+1], v_d_coulomb[i], v_d_coulomb[i+1], i);
            CSTable_CalcEnAttract[b][k].InitializeSpecificPoint(v_attract[i], v_attract[i+1], v_d_attract[i], v_d_attract[i+1], i);
            CSTable_CalcEnRepulse[b][k].InitializeSpecificPoint(v_repulse[i], v_repulse[i+1], v_d_repulse[i], v_d_repulse[i+1], i);
          }
          else {
            CSTable_CalcCoulomb[b][k].InitializeSpecificPoint(v_coulomb[i], v_coulomb[i], v_d_coulomb[i], v_d_coulomb[i], i);
            CSTable_CalcEnAttract[b][k].InitializeSpecificPoint(v_attract[i], v_attract[i], v_d_attract[i], v_d_attract[i], i);
            CSTable_CalcEnRepulse[b][k].InitializeSpecificPoint(v_repulse[i], v_repulse[i], v_d_repulse[i], v_d_repulse[i], i);
          }
        }
      }
    }
  }
}