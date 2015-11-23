#include "FFParticle.h"
#include "FFSetup.h" //For our setup info
#include "ConfigSetup.h"
#include "../lib/NumLib.h" //For Sq, Cb, and MeanA/G functions.

FFParticle::FFParticle() : mass(NULL), nameFirst(NULL), nameSec(NULL),
			   n(NULL), sigmaSq(NULL), epsilon_cn(NULL),
			   epsilon_cn_6(NULL), nOver6(NULL),
			   enCorrection(NULL), virCorrection(NULL),
			   shiftConst(NULL), An(NULL), Bn(NULL), Cn(NULL),
			   sig6(NULL), sign(NULL), rCut(0), rCutSq(0),
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
   sigmaSq = new double [size];
   isMartini = ffKind.isMARTINI;


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

   if(vdwKind == val.VDW_SHIFT_KIND)
     shiftConst = new double [size];
   
   rCut =  val.cutoff;
   rCutSq = rCut * rCut;
   if(vdwKind == val.VDW_SWITCH_KIND && isMartini)
   {
      An = new double [size];
      Bn = new double [size];
      Cn = new double [size];
      sign = new double [size];
      sig6 = new double [size];

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
	  double rRat = nbfix.sigma[i]/rCut, tc = 1.0;
	  num::Cb(sigmaSq[j], tc, nbfix.sigma[i]);
	  tc *= 0.5 * 4.0 * M_PI;
	  epsilon_cn[j] = nbfix.Cn[i] * nbfix.epsilon[i];
	  epsilon_cn_6[j] = epsilon_cn[j]*6;
	  nOver6[j] = n[j]/6;
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
#ifdef MIE_INT_ONLY
	    double repulse = num::POW(rRat2, rRat4, attract, n[j]);
#else
	    double repulse = pow(sqrt(rRat2), n[j]);
#endif
	    shiftConst[j] =  epsilon_cn[j] * (repulse-attract);
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

	  if(vdwKind == num::VDW_SHIFT_KIND)
	  {
	    double rRat2 = sigmaSq[idx]/rCutSq;
	    double rRat4 = rRat2 * rRat2;
	    double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
	    double repulse = num::POW(rRat2, rRat4, attract, n[j]);
#else
	    double repulse = pow(sqrt(rRat2), n[idx]);
#endif
	    shiftConst[idx] =  epsilon_cn[idx] * (repulse-attract);
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
	  }

      }
   }
}
