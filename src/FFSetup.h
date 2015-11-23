#ifndef FF_SETUP_H
#define FF_SETUP_H

#include <string> //for var names, etc.
#include <map> //for function handle storage.
#include <vector>
#include <sstream>

#include "InputAbstracts.h" //For readable base, etc.
#include "Reader.h" //For Reader object
#include "FFConst.h" //for forcefield constants
#include "../lib/BasicTypes.h" //for uint

namespace ff_setup
{
   extern const double KCAL_PER_MOL_TO_K; //503.21959899;
   extern const double RIJ_OVER_2_TO_SIG; //1.7817974362807;

   struct FFBase : SearchableBase
   {
      explicit FFBase(uint terms): SearchableBase(terms), numTerms(terms), 
	 multi(false) {}
      FFBase(uint terms, const bool mult): SearchableBase(terms), 
	 numTerms(terms), multi(mult) {}
      void IsCHARMM(const bool isCHRM) { isCHARMM = isCHRM; }
      std::string ReadKind(Reader & param, std::string const& firstKindName);
      
      double EnConvIfCHARMM(double val) const 
      { 
	 if (isCHARMM) 
	 { 
	    val *= KCAL_PER_MOL_TO_K;
	 }
	 return val;
      }
      std::string LoadLine(Reader & param, std::string const& firstVar);
      
      std::vector<std::string> name;
      uint numTerms;
      bool multi;
      bool isCHARMM;
   };

   struct Particle : ReadableBaseWithFirst, FFBase
   {
      std::vector<double> sigma, epsilon, sigma_1_4, epsilon_1_4;
      std::vector<uint> n, n_1_4;

      Particle(void) : FFBase(1) {}

      virtual void Read(Reader & param, std::string const& firstVar);
      void Add(double e, double s, const uint expN,
	       double e_1_4, double s_1_4, const uint expN_1_4);
#ifndef NDEBUG
      void PrintBrief();
#endif
   };

struct NBfix : ReadableBaseWithFirst, FFBase
   {
     std::vector<double> sigma, epsilon, Cn;
      std::vector<uint> n;

      NBfix() : FFBase(2) {}

      virtual void Read(Reader & param, std::string const& firstVar);
      void Add(double e, double s, const uint expN, double Cn);
   };


   struct Bond : ReadableBaseWithFirst, FFBase
   {
      std::vector<double> Kb, b0;
      //XXX This is not a real vector
      //XXX Do not use with std algorithms, they are not required to work
      std::vector<bool> fixed;

      Bond() : FFBase(2) {}
      virtual void Read(Reader & param, std::string const& firstVar);
      void Add(const double coeff, const double def);
      static const double FIXED;
#ifndef NDEBUG
      void PrintBrief();
#endif
   };

   struct Angle : ReadableBaseWithFirst, FFBase
   {
      std::vector<double> Ktheta, theta0, Kub, bUB0;
      //XXX This is not a real vector
      //XXX Do not use with std algorithms, they are not required to work
      std::vector<bool> hasUB;

      Angle() : FFBase(3) {}
      virtual void Read(Reader & param, std::string const& firstVar);
      void Add(const double coeff, const double def, const bool hsUB,
	       const double coeffUB, const double defUB);
#ifndef NDEBUG
      void PrintBrief();
#endif
   };

   struct Dihedral : ReadableBaseWithFirst, FFBase
   {
      std::map< std::string, std::vector<double> > Kchi, delta;
      std::map< std::string, std::vector<uint> > n;
      std::string last;
      uint countTerms;

      Dihedral(void): FFBase(4, true), last(""), countTerms(0) {}
      virtual void Read(Reader & param, std::string const& firstVar);
      void Add(std::string const& merged,
	       const double coeff, const uint index, const double def);
#ifndef NDEBUG
      void PrintBrief();
#endif
   };

   struct Improper : ReadableBaseWithFirst, FFBase
   {
      std::vector<double> Komega, omega0;
      
      Improper() : FFBase(4) {}
      virtual void Read(Reader & param, std::string const& firstVar);
      void Add(const double coeff, const double def);
#ifndef NDEBUG
      void PrintBrief();
#endif
   };
}

struct FFSetup
{
   ff_setup::Particle mie;
   ff_setup::NBfix nbfix;
   ff_setup::Bond bond;
   ff_setup::Angle angle;
   ff_setup::Dihedral dih;
   ff_setup::Improper imp;

   FFSetup(void) {}
   void Init(std::string const& fileName, const bool isCHARMM);

 private:  
   //Map variable names to functions 
   std::map<std::string, ReadableBaseWithFirst *> sectKind;  
   typedef std::map<std::string, ReadableBaseWithFirst *>::iterator sect_it;
   std::map<std::string, ReadableBaseWithFirst *>  SetReadFunctions
      (const bool isCHARMM);
   static const std::string paramFileAlias[];
   static const uint CHARMM_ALIAS_IDX;
   static const uint EXOTIC_ALIAS_IDX;
};


#endif /*FF_SETUP_H*/
