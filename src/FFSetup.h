/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_SETUP_H
#define FF_SETUP_H

#include <string> //for var names, etc.
#include <map> //for function handle storage.
#include <vector>
#include <sstream>

#include "InputAbstracts.h" //For readable base, etc.
#include "Reader.h" //For Reader object
#include "FFConst.h" //for forcefield constants
#include "BasicTypes.h" //for uint

namespace ff_setup
{
extern const real KCAL_PER_MOL_TO_K; //503.21959899;
extern const real RIJ_OVER_2_TO_SIG; //1.7817974362807;
extern const real RIJ_TO_SIG; //0.890898718

class FFBase : public SearchableBase
{
public:
  explicit FFBase(uint terms) : SearchableBase(terms), numTerms(terms),
    multi(false) {}
  FFBase(uint terms, const bool mult) : SearchableBase(terms),
    numTerms(terms), multi(mult) {}
  void setIsCHARMM(const bool isCHRM)
  {
    CHARMM = isCHRM;
  }
  bool isCHARMM(void) const
  {
    return CHARMM;
  }
  std::string ReadKind(Reader & param, std::string const& firstKindName);
  std::vector<std::string> & getnamelist()
  {
    return name;
  }

  bool validname(std::string & merged)
  {
    return (std::count(name.begin(), name.end(), merged) > 0);
  }
  std::string getname(const uint i) const
  {
    return name[i];
  }
  size_t getnamecnt() const
  {
    return name.size();
  }
  void clean_names()
  {
    std::vector<std::string>::iterator newEnd;
    newEnd = std::unique(name.begin(), name.end());
    name.erase(newEnd, name.end());
  }

  real EnConvIfCHARMM(real val) const
  {
    if (CHARMM) {
      val *= KCAL_PER_MOL_TO_K;
    }
    return val;
  }

  std::string LoadLine(Reader & param, std::string const& firstVar);

//	private:
  std::vector<std::string> name;
  uint numTerms;
  bool multi;
private:
  bool CHARMM;
};

class Particle : public ReadableBaseWithFirst, public FFBase
{
public:
  Particle(void) : FFBase(1) {}

  virtual void Read(Reader & param, std::string const& firstVar);
  void Add(real e, real s, const uint expN,
           real e_1_4, real s_1_4, const uint expN_1_4);
#ifndef NDEBUG
  void PrintBrief();
#endif
//	private:
  std::vector<real> sigma, epsilon, sigma_1_4, epsilon_1_4;
  std::vector<uint> n, n_1_4;

};

class NBfix : public ReadableBaseWithFirst, public FFBase
{
public:
  NBfix() : FFBase(2) {}

  virtual void Read(Reader & param, std::string const& firstVar);
  void Add(real e, real s,
#ifdef MIE_INT_ONLY
           const uint expN,
#else
           const real expN,
#endif
           real e_1_4, real s_1_4,
#ifdef MIE_INT_ONLY
           const uint expN_1_4
#else
           const real expN_1_4
#endif
          );
//	private:
  std::vector<real> sigma, epsilon, sigma_1_4, epsilon_1_4;
#ifdef MIE_INT_ONLY
  std::vector<uint> n, n_1_4;
#else
  std::vector<real> n, n_1_4;
#endif
};


class Bond : public ReadableBaseWithFirst, public FFBase
{
public:
  Bond() : FFBase(2) {}
  virtual void Read(Reader & param, std::string const& firstVar);
  void Add(const real coeff, const real def);
  real *CopyKb() const
  {
    return vect::transfer(Kb);
  }
  real *Copyb0() const
  {
    return vect::transfer(b0);
  }
  bool * Copyfixed() const
  {
    return vect::transfer(fixed);
  }
  size_t getKbcnt() const
  {
    return Kb.size();
  }
  real GetKb(uint i) const
  {
    return Kb[i];
  }
  real Getb0(uint i) const
  {
    return b0[i];
  }
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  static const real FIXED;
  std::vector<real> Kb, b0;
  //XXX This is not a real vector
  //XXX Do not use with std algorithms, they are not required to work
  std::vector<bool> fixed;
};

class Angle : public ReadableBaseWithFirst, public FFBase
{
public:
  Angle() : FFBase(3) {}
  virtual void Read(Reader & param, std::string const& firstVar);
  void Add(const real coeff, const real def, const bool hsUB,
           const real coeffUB, const real defUB);
  real *CopyKtheta() const
  {
    return vect::transfer(Ktheta);
  }
  bool * Copyfixed() const
  {
    return vect::transfer(fixed);
  }
  real *Copytheta0() const
  {
    return vect::transfer(theta0);
  }
  size_t getKthetacnt() const
  {
    return Ktheta.size();
  }
  real GetKtheta(uint i) const
  {
    return Ktheta[i];
  }
  real Gettheta0(uint i) const
  {
    return theta0[i];
  }
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  static const real FIXED;
  std::vector<real> Ktheta, theta0, Kub, bUB0;
  //XXX This is not a real vector
  //XXX Do not use with std algorithms, they are not required to work
  std::vector<bool> hasUB;
  std::vector<bool> fixed;
};

class Dihedral : public ReadableBaseWithFirst, public FFBase
{
public:
  Dihedral(void) : FFBase(4, true), last(""), countTerms(0) {}
  virtual void Read(Reader & param, std::string const& firstVar);
  void Add(std::string const& merged,
           const real coeff, const uint index, const real def);
  uint getTerms() const
  {
    return countTerms;
  }
  uint append(std::string & s, real * Kchi_in, real * delta_in, uint * n_in, uint count) const
  {
    std::map<std::string, std::vector<uint> >::const_iterator itUInt = n.find(s);
    std::copy(itUInt->second.begin(), itUInt->second.end(), n_in + count);
    std::map<std::string, std::vector<real> >::const_iterator itDbl = Kchi.find(s);
    std::copy(itDbl->second.begin(), itDbl->second.end(), Kchi_in + count);
    itDbl = delta.find(s);
    std::copy(itDbl->second.begin(), itDbl->second.end(), delta_in + count);
    return itDbl->second.size();
  }

  uint GetSizeDih(std::string s) const
  {
    std::map< std::string, std::vector<real> >::const_iterator it = Kchi.find(s);
    return it->second.size();
  }

  real GetKchi(std::string s, uint pos) const
  {
    std::map< std::string, std::vector<real> >::const_iterator it = Kchi.find(s);
    return it->second[pos];
  }

  real Getdelta(std::string s, uint pos) const
  {
    std::map< std::string, std::vector<real> >::const_iterator it = delta.find(s);
    return it->second[pos];
  }

  uint Getn(std::string s, uint pos) const
  {
    std::map< std::string, std::vector<uint> >::const_iterator it = n.find(s);
    return it->second[pos];
  }

#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  std::map< std::string, std::vector<real> > Kchi, delta;
  std::map< std::string, std::vector<uint> > n;
  std::string last;
  uint countTerms;
};

class Improper : public ReadableBaseWithFirst, public FFBase
{
public:
  Improper() : FFBase(4) {}
  virtual void Read(Reader & param, std::string const& firstVar);
  void Add(const real coeff, const real def);
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  std::vector<real> Komega, omega0;
};
}

class FFSetup
{
public:
  FFSetup(void) {}
  void Init(std::string const& fileName, const bool isCHARMM);

  ff_setup::Particle mie;
  ff_setup::NBfix nbfix;
  ff_setup::Bond bond;
  ff_setup::Angle angle;
  ff_setup::Dihedral dih;
  ff_setup::Improper imp;

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
