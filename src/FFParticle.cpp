/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFParticle.h"
#include "FFSetup.h" //For our setup info
#include "ConfigSetup.h"
#include "../lib/NumLib.h" //For Sq, Cb, and MeanA/G functions.

FFParticle::FFParticle() : mass(NULL), nameFirst(NULL), nameSec(NULL),
			   n(NULL), n_1_4(NULL), sigmaSq(NULL),
			   sigmaSq_1_4(NULL),epsilon_cn(NULL),
			   epsilon_cn_1_4(NULL), epsilon_cn_6(NULL),
			   epsilon_cn_6_1_4(NULL), nOver6(NULL),
			   nOver6_1_4(NULL), enCorrection(NULL),
			   virCorrection(NULL), shiftConst(NULL), An(NULL),
			   Bn(NULL), Cn(NULL), shiftConst_1_4(NULL), An_1_4(NULL), Bn_1_4(NULL),
			   Cn_1_4(NULL), sig6(NULL), sign(NULL),
			   sig6_1_4(NULL), sign_1_4(NULL), rCut(0), rCutSq(0),
			   rOnSq(0), rOn(0), A6(0), B6(0), C6(0), factor1(0),
			   factor2(0) {}

FFParticle::~FFParticle(void)
{
   delete[] mass;
   delete[] nameFirst;
   delete[] nameSec;

   delete[] sigmaSq;
   delete[] n;
   delete[] epsilon_cn;
   delete[] epsilon_cn_6;
   delete[] nOver6;
   // parameter for 1-4 interaction, same one will be used for 1-3 interaction
   delete[] sigmaSq_1_4;
   delete[] n_1_4;
   delete[] epsilon_cn_1_4;
   delete[] epsilon_cn_6_1_4;
   delete[] nOver6_1_4;

   delete[] enCorrection;
   delete[] virCorrection;

   delete[] shiftConst;
   delete[] An;
   delete[] Bn;
   delete[] Cn;
   delete[] sig6;
   delete[] sign;
   // parameter for 1-4 interaction, same one will be used for 1-3 interaction
   delete[] shiftConst_1_4;
   delete[] An_1_4;
   delete[] Bn_1_4;
   delete[] Cn_1_4;
   delete[] sig6_1_4;
   delete[] sign_1_4;
}

void FFParticle::Init(ff_setup::Particle const& mie,
		      ff_setup::NBfix const& nbfix,
		      config_setup::FFValues const& val,
		      config_setup::FFKind const& ffKind)
{
   count = mie.epsilon.size(); //Get # particles read
   //Size LJ particle kind arrays
   mass = new double [count];
   vdwKind = val.VDW_KIND;
   
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
   epsilon_cn = new double [size];
   epsilon_cn_6 = new double [size];   
   nOver6 = new double [size]; 
   sigmaSq = new double [size];

   epsilon_cn_1_4 = new double [size];
   epsilon_cn_6_1_4 = new double [size];
   nOver6_1_4 = new double [size]; 
   sigmaSq_1_4 = new double [size];

   enCorrection = new double [size];
   virCorrection = new double [size];

   if(vdwKind == val.VDW_SHIFT_KIND)
   {
     shiftConst = new double [size];
     shiftConst_1_4 = new double [size];
   }
   
   rCut =  val.cutoff;
   rCutSq = rCut * rCut;
   if(vdwKind == val.VDW_SWITCH_KIND && isMartini)
   {
      An = new double [size];
      Bn = new double [size];
      Cn = new double [size];
      An_1_4 = new double [size];
      Bn_1_4 = new double [size];
      Cn_1_4 = new double [size];
      sign = new double [size];
      sig6 = new double [size];
      sign_1_4 = new double [size];
      sig6_1_4 = new double [size];

      rOn = val.rswitch;
      rOnSq = rOn * rOn;
      A6 = 6.0 * ((6.0+1)*rOn-(6.0+4)*rCut)/(pow(rCut,6.0+2)*
						 pow(rCut-rOn, 2));
      B6 = -6.0 * ((6.0+1)*rOn-(6.0+3)*rCut)/(pow(rCut,6.0+2)*
						 pow(rCut-rOn, 3));
      C6 = 1.0/pow(rCut, 6.0)-A6/3.0*pow(rCut-rOn,3)-B6/4.0*
	pow(rCut-rOn,4);
   }
   if(vdwKind == val.VDW_SWITCH_KIND && !isMartini)
   {
     rOn = val.rswitch;
     rOnSq = rOn * rOn;
     
     factor1 = rCutSq - 3 * rOnSq;
     factor2 = pow((rCutSq - rOnSq), -3);
   }
   
   Blend(mie, rCut);
   AdjNBfix(mie, nbfix, rCut);
}

double FFParticle::EnergyLRC(const uint kind1, const uint kind2) const
{ return enCorrection[FlatIndex(kind1, kind2)]; }


double FFParticle::VirialLRC(const uint kind1, const uint kind2) const
{ return virCorrection[FlatIndex(kind1, kind2)]; }

void FFParticle::AdjNBfix(ff_setup::Particle const& mie,
		 ff_setup::NBfix const& nbfix, const double rCut)
{
   uint size = num::Sq(count);
   for(uint i = 0; i < nbfix.epsilon.size(); i++)
   {
     for(uint j = 0; j < size; j++)
     {
       if(nbfix.name[i] == nameFirst[j] || nbfix.name[i] ==  nameSec[j])
       {
	  n[j] = nbfix.n[i];
	  n_1_4[j] = nbfix.n_1_4[i];
	  double rRat = nbfix.sigma[i]/rCut, tc = 1.0;
	  num::Cb(sigmaSq[j], tc, nbfix.sigma[i]);
	  num::Cb(sigmaSq_1_4[j], tc, nbfix.sigma_1_4[i]);
	  tc *= 0.5 * 4.0 * M_PI;
	  double cn = n[j]/(n[j]-6) *pow(n[j]/6, (6/(n[j]-6)));
	  double cn_1_4 = n_1_4[j]/(n_1_4[j]-6) *
	   pow(n_1_4[j]/6, (6/(n_1_4[j]-6)));
	  epsilon_cn[j] = cn * nbfix.epsilon[i];
	  epsilon_cn_1_4[j] = cn_1_4 * nbfix.epsilon_1_4[i];
	  epsilon_cn_6[j] = epsilon_cn[j]*6;
	  epsilon_cn_6_1_4[j] = epsilon_cn_1_4[j]*6;
	  nOver6[j] = n[j]/6;
	  nOver6_1_4[j] = n_1_4[j]/6;
	  enCorrection[j] = tc/(n[j]-3) * epsilon_cn[j] * 
            ( pow(rRat, n[j]-3) - 
	      (double)(n[j]-3.0)/3.0 * pow(rRat, 3) );
          virCorrection[j] = tc/(n[j]-3) * epsilon_cn_6[j] *
            ( (double)(n[j])/6.0 * pow(rRat, n[j]-3) - 
	      (double)(n[j]-3.0)/3.0 * pow(rRat, 3) );


	  if(vdwKind == num::VDW_SHIFT_KIND)
	  {
	    double rRat2 = sigmaSq[j]/rCutSq;
	    double rRat4 = rRat2 * rRat2;
	    double attract = rRat4 * rRat2;
	    //for 1-4 interaction
	    double rRat2_1_4 = sigmaSq_1_4[j]/rCutSq;
	    double rRat4_1_4 = rRat2_1_4 * rRat2_1_4;
	    double attract_1_4 = rRat4_1_4 * rRat2_1_4;
#ifdef MIE_INT_ONLY
	    double repulse = num::POW(rRat2, rRat4, attract, n[j]);
	    double repulse_1_4 =
	      num::POW(rRat2_1_4, rRat4_1_4, attract_1_4, n_1_4[j]);
#else
	    double repulse = pow(sqrt(rRat2), n[j]);
	    double repulse_1_4 = pow(sqrt(rRat2_1_4), n_1_4[j]);
#endif
	    shiftConst[j] =  epsilon_cn[j] * (repulse - attract);
	    shiftConst_1_4[j] =  epsilon_cn_1_4[j] *
	      (repulse_1_4 - attract_1_4);
	  }

	  if(vdwKind == num::VDW_SWITCH_KIND && isMartini)
	  {
	    double pn = n[j];
	    An[j] = pn * ((pn+1)*rOn-(pn+4)*rCut)/(pow(rCut,pn+2)*
						     pow(rCut-rOn,2));
	    Bn[j] = -pn * ((pn+1)*rOn-(pn+3)*rCut)/(pow(rCut,pn+2)*
						      pow(rCut-rOn,3));
	    Cn[j] = 1.0/pow(rCut,pn)-An[j]/3.0*pow(rCut-rOn,3)-
	      Bn[j]/4.0*pow(rCut-rOn,4);
	    sig6[j] = pow(nbfix.sigma[i], 6);
	    sign[j] = pow(nbfix.sigma[i], pn);

	    // for 1-4 interaction
	    double pn_1_4 = n_1_4[j];
	    An_1_4[j] = pn_1_4 * ((pn_1_4 + 1)*rOn-(pn_1_4 + 4)*rCut)/
	      (pow(rCut,pn_1_4 + 2)*pow(rCut-rOn,2));
	    Bn_1_4[j] = -pn_1_4 * ((pn_1_4 + 1)*rOn-(pn_1_4 + 3)*rCut)/
	      (pow(rCut,pn_1_4+2)*pow(rCut-rOn,3));
	    Cn_1_4[j] = 1.0/pow(rCut,pn_1_4)-An_1_4[j]/3.0*pow(rCut-rOn,3)-
	      Bn_1_4[j]/4.0*pow(rCut-rOn,4);
	    sig6_1_4[j] = pow(nbfix.sigma_1_4[i], 6);
	    sign_1_4[j] = pow(nbfix.sigma_1_4[i], pn_1_4);
	  }
       }
     }
   }    
}

void FFParticle::Blend(ff_setup::Particle const& mie,
		       const double rCut)
{
   for(uint i = 0; i < count; ++i)
   {

      for(uint j = 0; j < count; ++j)
      {
	 uint idx = FlatIndex(i, j);
	 //get all name combination for using in nbfix
	 nameFirst[idx] = mie.name[i]; 
	 nameFirst[idx] += mie.name[j];
	 nameSec[idx] = mie.name[j]; 
	 nameSec[idx] += mie.name[i];

	 n[idx] = num::MeanA(mie.n, mie.n, i, j);
	 n_1_4[idx] = num::MeanA(mie.n_1_4, mie.n_1_4, i, j);
	 double cn = n[idx]/(n[idx]-6) * pow(n[idx]/6, (6/(n[idx]-6))); 
	 double sigma = num::MeanA(mie.sigma, mie.sigma, i, j);
	 double cn_1_4 = n_1_4[idx]/(n_1_4[idx]-6) *
	   pow(n_1_4[idx]/6, (6/(n_1_4[idx]-6))); 
	 double sigma_1_4 = num::MeanA(mie.sigma_1_4, mie.sigma_1_4, i, j);
	 double tc = 1.0;
	 double rRat = sigma/rCut; 

	 num::Cb(sigmaSq[idx], tc, sigma);
	 num::Cb(sigmaSq_1_4[idx], tc, sigma_1_4);
	 tc *= 0.5 * 4.0 * M_PI;
	 epsilon_cn[idx] = 
	    cn * num::MeanG(mie.epsilon, mie.epsilon, i, j);
	 epsilon_cn_1_4[idx] = 
	    cn * num::MeanG(mie.epsilon_1_4, mie.epsilon_1_4, i, j);
	 epsilon_cn_6[idx] = epsilon_cn[idx]*6;
	 epsilon_cn_6_1_4[idx] = epsilon_cn_1_4[idx]*6;
	 nOver6[idx] = n[idx]/6;
	 nOver6_1_4[idx] = n_1_4[idx]/6;
	 enCorrection[idx] = tc/(n[idx]-3) * epsilon_cn[idx] * 
            ( pow(rRat, n[idx]-3) - 
	      (double)(n[idx]-3.0)/3.0 * pow(rRat, 3) );
         virCorrection[idx] = tc/(n[idx]-3) * epsilon_cn_6[idx] *
            ( (double)(n[idx])/6.0 * pow(rRat, n[idx]-3) - 
	      (double)(n[idx]-3.0)/3.0 * pow(rRat, 3) );

	  if(vdwKind == num::VDW_SHIFT_KIND)
	  {
	    double rRat2 = sigmaSq[idx]/rCutSq;
	    double rRat4 = rRat2 * rRat2;
	    double attract = rRat4 * rRat2;
	    //for 1-4 interaction
	    double rRat2_1_4 = sigmaSq_1_4[idx]/rCutSq;
	    double rRat4_1_4 = rRat2_1_4 * rRat2_1_4;
	    double attract_1_4 = rRat4_1_4 * rRat2_1_4;
	    
#ifdef MIE_INT_ONLY
	    double repulse = num::POW(rRat2, rRat4, attract, n[idx]);
	    double repulse_1_4 =
	      num::POW(rRat2_1_4, rRat4_1_4, attract_1_4, n_1_4[idx]);
#else
	    double repulse = pow(sqrt(rRat2), n[idx]);
	    double repulse_1_4 = pow(sqrt(rRat2_1_4), n_1_4[idx]);
#endif
	    shiftConst[idx] =  epsilon_cn[idx] * (repulse - attract);
	    shiftConst_1_4[idx] =  epsilon_cn_1_4[idx] *
	      (repulse_1_4 - attract_1_4);
	  }
	  
	  if(vdwKind == num::VDW_SWITCH_KIND && isMartini)
	  {
	    double pn = n[idx];
	    An[idx] = pn * ((pn+1)*rOn-(pn+4)*rCut)/(pow(rCut,pn+2)*
						     pow(rCut-rOn,2));
	    Bn[idx] = -pn * ((pn+1)*rOn-(pn+3)*rCut)/(pow(rCut,pn+2)*
						      pow(rCut-rOn,3));
	    Cn[idx] = 1.0/pow(rCut,pn)-An[idx]/3.0*pow(rCut-rOn,3)-
	      Bn[idx]/4.0*pow(rCut-rOn,4);
	    sig6[idx] = pow(sigma, 6);
	    sign[idx] = pow(sigma, pn);

	    // for 1-4 interaction
	    double pn_1_4 = n_1_4[idx];
	    An_1_4[idx] = pn_1_4 * ((pn_1_4+1)*rOn-(pn_1_4+4)*rCut)/
	      (pow(rCut,pn_1_4+2)*pow(rCut-rOn,2));
	    Bn_1_4[idx] = -pn_1_4 * ((pn_1_4+1)*rOn-(pn_1_4+3)*rCut)/
	      (pow(rCut,pn_1_4+2)*pow(rCut-rOn,3));
	    Cn_1_4[idx] = 1.0/pow(rCut,pn_1_4)-An_1_4[idx]/3.0*pow(rCut-rOn,3)-
	      Bn_1_4[idx]/4.0*pow(rCut-rOn,4);
	    sig6_1_4[idx] = pow(sigma_1_4, 6);
	    sign_1_4[idx] = pow(sigma_1_4, pn_1_4);
	  }

      }
   }
}

