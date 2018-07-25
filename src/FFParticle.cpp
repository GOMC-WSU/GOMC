/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFParticle.h"
#include "ConfigSetup.h"
#include "NumLib.h" //For Sq, Cb, and MeanA/G functions.
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

FFParticle::FFParticle() : mass(NULL), nameFirst(NULL), nameSec(NULL),
  n(NULL), n_1_4(NULL), sigmaSq(NULL), sigmaSq_1_4(NULL), epsilon_cn(NULL),
  epsilon(NULL), epsilon_1_4(NULL), epsilon_cn_1_4(NULL), epsilon_cn_6(NULL),
  epsilon_cn_6_1_4(NULL), nOver6(NULL), nOver6_1_4(NULL), enCorrection(NULL),
  virCorrection(NULL) 
#ifdef GOMC_CUDA
  , varCUDA(NULL)
#endif 
  {
    rCut = rCutSq = 0.0;
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
                      ff_setup::NBfix const& nbfix,
                      config_setup::SystemVals const& sys,
                      config_setup::FFKind const& ffKind)
{
#ifdef GOMC_CUDA
  // Variables for GPU stored in here
  varCUDA = new VariablesCUDA();
#endif
  count = mie.epsilon.size(); //Get # particles read
  //Size LJ particle kind arrays
  mass = new double [count];
  vdwKind = sys.ff.VDW_KIND;
  vdwGeometricSigma = sys.ff.vdwGeometricSigma;

  //Size LJ-LJ pair arrays
  uint size = num::Sq(count);
  nameFirst = new std::string [size];
  nameSec = new std::string [size];
  isMartini = ffKind.isMARTINI;


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

  rCut =  sys.ff.cutoff;
  rCutSq = rCut * rCut;
  rCutLow = sys.ff.cutoffLow;
  rCutLowSq = rCutLow * rCutLow;
  scaling_14 = sys.elect.oneFourScale;
  ewald = sys.elect.ewald;
  alpha = sys.elect.alpha;

  //Combining VDW parameter
  Blend(mie, rCut);
  //Adjusting VDW parameter using NBFIX
  AdjNBfix(mie, nbfix, rCut);

#ifdef GOMC_CUDA
  double diElectric_1 = 1.0 / sys.elect.dielectric;
  double rOn = sys.ff.rswitch;
  InitGPUForceField(*varCUDA, sigmaSq, epsilon_cn, n, vdwKind, isMartini,
                    count, rCut, rCutLow, rOn, alpha, ewald, diElectric_1);
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
      
      if(vdwGeometricSigma) {
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
