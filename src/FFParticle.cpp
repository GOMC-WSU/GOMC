/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFParticle.h"
#include "FFSetup.h" //For our setup info
#include "../lib/NumLib.h" //For Sq, Cb, and MeanA/G functions.

FFParticle::FFParticle() : mass(NULL), name(NULL), n(NULL), sigmaSq(NULL), 
   epsilon_cn(NULL), epsilon_cn_6(NULL), nOver6(NULL), enCorrection(NULL),
   virCorrection(NULL) {}

FFParticle::~FFParticle(void)
{
   /*delete[] mass;
   delete[] name;
   delete[] sigmaSq;
   delete[] n;
   delete[] epsilon_cn;
   delete[] epsilon_cn_6;
   delete[] nOver6;
   delete[] sigmaSq_1_4;
   delete[] n_1_4;
   delete[] epsilon_cn_1_4;
   delete[] epsilon_cn_6_1_4;
   delete[] nOver6_1_4;
   delete[] enCorrection;
   delete[] virCorrection;*/
}

void FFParticle::Init(ff_setup::Particle const& mie, const double rCut)
{
   count = mie.epsilon.size(); //Get # particles read
   //Size LJ particle kind arrays
   mass = new double [count];
   name = new char [count*ff::part::nm_len];
   //Size LJ-LJ pair arrays
   uint size = num::Sq(count);
   sigmaSq = new double [size];
#ifdef MIE_INT_ONLY
   n = new uint [size];
#else 
   n = new double [size];
#endif
   epsilon_cn = new double [size];
   epsilon_cn_6 = new double [size];
   nOver6 = new double [size]; 
   enCorrection = new double [size];
   virCorrection = new double [size];
   
   Blend(mie, rCut);
}


double FFParticle::EnergyLRC(const uint kind1, const uint kind2) const
{ return enCorrection[FlatIndex(kind1, kind2)]; }


double FFParticle::VirialLRC(const uint kind1, const uint kind2) const
{ return virCorrection[FlatIndex(kind1, kind2)]; }

void FFParticle::Blend(ff_setup::Particle const& mie,
		       const double rCut)
{
   for(uint i = 0; i < count; ++i)
   {

      for(uint j = 0; j < count; ++j)
      {
	 uint idx = FlatIndex(i, j);
	 n[idx] = num::MeanA(mie.n, mie.n, i, j);
	 double cn = n[idx]/(n[idx]-6) * pow(n[idx]/6, (6/(n[idx]-6))), 
	    sigma = num::MeanA(mie.sigma, mie.sigma, i, j), tc = 1.0,
	    rRat = sigma/rCut; 
	 num::Cb(sigmaSq[idx], tc, sigma);
	 tc *= 0.5 * 4.0 * M_PI;
	 epsilon_cn[idx] = 
	    cn * num::MeanG(mie.epsilon, mie.epsilon, i, j);
	 epsilon_cn_6[idx] = epsilon_cn[idx]*6;
	 nOver6[idx] = n[idx]/6;
	 enCorrection[idx] = tc/(n[idx]-3) * epsilon_cn[idx] * 
            ( pow(rRat, n[idx]-3) - 
	      (double)(n[idx]-3.0)/3.0 * pow(rRat, 3) );
         virCorrection[idx] = tc/(n[idx]-3) * epsilon_cn_6[idx] *
            ( (double)(n[idx])/6.0 * pow(rRat, n[idx]-3) - 
	      (double)(n[idx]-3.0)/3.0 * pow(rRat, 3) );
    }
  }
}

