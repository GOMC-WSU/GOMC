/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FF_SETUP_H
#define FF_SETUP_H

#include <map> //for function handle storage.
#include <sstream>
#include <string> //for var names, etc.
#include <vector>

#include "BasicTypes.h"     //for uint
#include "ConfigSetup.h"    //for access to structure
#include "FFConst.h"        //for forcefield constants
#include "InputAbstracts.h" //For readable base, etc.
#include "InputFileReader.h"
#include "Reader.h" //For Reader object

namespace config_setup {
struct FileName;
}

namespace ff_setup {
extern const double KCAL_PER_MOL_TO_K; // 503.21959899;
extern const double RIJ_OVER_2_TO_SIG; // 1.7817974362807;
extern const double RIJ_TO_SIG;        // 0.890898718

class FFBase : public SearchableBase {
public:
  explicit FFBase(uint terms)
      : SearchableBase(terms), numTerms(terms), multi(false) {}
  FFBase(uint terms, const bool mult)
      : SearchableBase(terms), numTerms(terms), multi(mult) {}
  void setIsCHARMM(const bool isCHRM) { CHARMM = isCHRM; }
  bool isCHARMM(void) const { return CHARMM; }
  std::string ReadKind(Reader &param, std::string const &firstKindName);
  std::vector<std::string> &getnamelist() { return name; }

  bool validname(const std::string &merged) const {
    return (std::count(name.begin(), name.end(), merged) > 0);
  }
  std::string getname(const uint i) const { return name[i]; }
  size_t getnamecnt() const { return name.size(); }
  void clean_names() {
    std::vector<std::string>::iterator newEnd;
    newEnd = std::unique(name.begin(), name.end());
    name.erase(newEnd, name.end());
  }

  double EnConvIfCHARMM(double val) const {
    if (CHARMM) {
      val *= KCAL_PER_MOL_TO_K;
    }
    return val;
  }

  std::string LoadLine(Reader &param, std::string const &firstVar);

  //    private:
  std::vector<std::string> name;
  uint numTerms;
  bool multi;

private:
  bool CHARMM;
};

class Particle : public ReadableBaseWithFirst, public FFBase {
public:
  Particle(void) : FFBase(1) {}

  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(double e, double s, const double expN, double e_1_4, double s_1_4,
           const double expN_1_4);
#ifndef NDEBUG
  void PrintBrief();
#endif
  //    private:
  std::vector<double> sigma, epsilon, sigma_1_4, epsilon_1_4;
  std::vector<double> n, n_1_4;
};

class NBfix : public ReadableBaseWithFirst, public FFBase {
public:
  NBfix() : FFBase(2) {}

  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(double e, double s, const double expN, double e_1_4, double s_1_4,
           const double expN_1_4);
  //    private:
  std::vector<double> sigma, epsilon, sigma_1_4, epsilon_1_4;
  std::vector<double> n, n_1_4;
};

class Bond : public ReadableBaseWithFirst, public FFBase {
public:
  Bond() : FFBase(2) {}
  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(const double coeff, const double def);
  double *CopyKb() const { return vect::transfer(Kb); }
  double *Copyb0() const { return vect::transfer(b0); }
  bool *Copyfixed() const { return vect::transfer(fixed); }
  size_t getKbcnt() const { return Kb.size(); }
  double GetKb(uint i) const { return Kb[i]; }
  double Getb0(uint i) const { return b0[i]; }
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  static const double FIXED;
  std::vector<double> Kb, b0;
  // XXX This is not a real vector
  // XXX Do not use with std algorithms, they are not required to work
  std::vector<bool> fixed;
};

class Angle : public ReadableBaseWithFirst, public FFBase {
public:
  Angle() : FFBase(3) {}
  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(const double coeff, const double def, const bool hsUB,
           const double coeffUB, const double defUB);
  double *CopyKtheta() const { return vect::transfer(Ktheta); }
  bool *Copyfixed() const { return vect::transfer(fixed); }
  double *Copytheta0() const { return vect::transfer(theta0); }
  size_t getKthetacnt() const { return Ktheta.size(); }
  double GetKtheta(uint i) const { return Ktheta[i]; }
  double Gettheta0(uint i) const { return theta0[i]; }
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  static const double FIXED;
  std::vector<double> Ktheta, theta0, Kub, bUB0;
  // XXX This is not a real vector
  // XXX Do not use with std algorithms, they are not required to work
  std::vector<bool> hasUB;
  std::vector<bool> fixed;
};

class Dihedral : public ReadableBaseWithFirst, public FFBase {
public:
  Dihedral(void) : FFBase(4, true), last(""), countTerms(0) {}
  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(std::string const &merged, const double coeff, const uint index,
           const double def);
  uint getTerms() const { return countTerms; }
  uint append(std::string &s, double *Kchi_in, double *delta_in, uint *n_in,
              uint count) const {
    std::map<std::string, std::vector<uint>>::const_iterator itUInt = n.find(s);
    std::copy(itUInt->second.begin(), itUInt->second.end(), n_in + count);
    std::map<std::string, std::vector<double>>::const_iterator itDbl =
        Kchi.find(s);
    std::copy(itDbl->second.begin(), itDbl->second.end(), Kchi_in + count);
    itDbl = delta.find(s);
    std::copy(itDbl->second.begin(), itDbl->second.end(), delta_in + count);
    return itDbl->second.size();
  }

  uint GetSizeDih(std::string s) const {
    std::map<std::string, std::vector<double>>::const_iterator it =
        Kchi.find(s);
    return it->second.size();
  }

  double GetKchi(std::string s, uint pos) const {
    std::map<std::string, std::vector<double>>::const_iterator it =
        Kchi.find(s);
    return it->second[pos];
  }

  double Getdelta(std::string s, uint pos) const {
    std::map<std::string, std::vector<double>>::const_iterator it =
        delta.find(s);
    return it->second[pos];
  }

  uint Getn(std::string s, uint pos) const {
    std::map<std::string, std::vector<uint>>::const_iterator it = n.find(s);
    return it->second[pos];
  }

#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  std::map<std::string, std::vector<double>> Kchi, delta;
  std::map<std::string, std::vector<uint>> n;
  std::string last;
  uint countTerms;
};

class Improper : public ReadableBaseWithFirst, public FFBase {
public:
  Improper() : FFBase(4) {}
  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(const double coeff, const double def);
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  std::vector<double> Komega, omega0;
};

class CMap : public ReadableBaseWithFirst, public FFBase {
public:
  CMap() : FFBase(4) {}
  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(const double coeff, const double def);
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  std::vector<double> Komega, omega0;
};

class HBond : public ReadableBaseWithFirst, public FFBase {
public:
  HBond() : FFBase(4) {}
  virtual void Read(Reader &param, std::string const &firstVar);
  void Add(const double coeff, const double def);
#ifndef NDEBUG
  void PrintBrief();
#endif
private:
  std::vector<double> Komega, omega0;
};
} // namespace ff_setup

class FFSetup {
public:
  FFSetup(void) {}
  void Init(const std::vector<config_setup::FileName> &fileName,
            const bool isCHARMM);

  ff_setup::Particle mie;
  ff_setup::NBfix nbfix;
  ff_setup::Bond bond;
  ff_setup::Angle angle;
  ff_setup::Dihedral dih;
  ff_setup::Improper imp;
  ff_setup::CMap cmap;
  ff_setup::HBond hbond;

private:
  // Map variable names to functions
  std::map<std::string, ReadableBaseWithFirst *> sectKind;
  typedef std::map<std::string, ReadableBaseWithFirst *>::iterator sect_it;
  std::map<std::string, ReadableBaseWithFirst *>
  SetReadFunctions(const bool isCHARMM);
  static const std::string paramFileAlias[];
  static const uint CHARMM_ALIAS_IDX;
  static const uint EXOTIC_ALIAS_IDX;
  bool hasEnding(std::string const &fullString, std::string const &ending);
};

#endif /*FF_SETUP_H*/
