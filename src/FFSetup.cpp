/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFSetup.h"
#include <algorithm>
#include <iostream>

#ifndef NDEBUG
#include <numeric>
#endif

#include "GeomLib.h"


const uint FFSetup::CHARMM_ALIAS_IDX = 0;
const uint FFSetup::EXOTIC_ALIAS_IDX = 1;
const std::string FFSetup::paramFileAlias[] =
{"CHARMM-Style parameter file", "EXOTIC-Style parameter file"};
const real ff_setup::KCAL_PER_MOL_TO_K = 503.21959899;
const real ff_setup::RIJ_OVER_2_TO_SIG = 1.7817974362807;
const real ff_setup::RIJ_TO_SIG = 0.890898718;

const real ff_setup::Bond::FIXED = 99999999;
const real ff_setup::Angle::FIXED = 99999999;

//Map variable names to functions
std::map<std::string, ReadableBaseWithFirst *>
FFSetup::SetReadFunctions(const bool isCHARMM)
{
  std::map<std::string, ReadableBaseWithFirst *> funct;
  //From CHARMM style file.
  funct["BONDS"] = &bond;
  funct["ANGLES"] = &angle;
  funct["DIHEDRALS"] = &dih;
  funct["IMPROPER"] = &imp;
  if (isCHARMM) {
    funct["NONBONDED"] = &mie;
    funct["NBFIX"] = &nbfix;
  } else {
    //Unique to exotic file
    funct["NONBONDED_MIE"] = &mie;
    funct["NBFIX_MIE"] = &nbfix;
  }
  for (sect_it it = funct.begin(); it != funct.end(); ++it) {
    (dynamic_cast<ff_setup::FFBase *>(it->second))->setIsCHARMM(isCHARMM);
  }
  return funct;
}

void FFSetup::Init(std::string const& name, const bool isCHARMM)
{
  using namespace std;
  sectKind = SetReadFunctions(isCHARMM);
  string currSectName = "", varName = "";
  string commentChar = "*!";
  string commentStr = "REMARK set AEXP REXP HAEX AAEX NBOND "
                      "CUTNB END CTONN EPS VSWI NBXM INHI";
  map<string, ReadableBaseWithFirst *>::const_iterator sect, currSect;

  Reader param(name,
               paramFileAlias[isCHARMM ? CHARMM_ALIAS_IDX : EXOTIC_ALIAS_IDX],
               true, &commentStr, true, &commentChar);
  param.open();
  while (param.Read(varName)) {
    sect = sectKind.find(varName);
    if ( sect != sectKind.end() ) {
      param.SkipLine(); //Skip rest of line for sect. heading
      currSectName = varName;
      currSect = sect; //Save for later calls.
      std::cout << "Reading " << currSectName << " parameters.\n";
    } else
      currSect->second->Read(param, varName);
  }

  param.close();

  //check if we read nonbonded parameter
  if(mie.sigma.size() == 0) {
    if(isCHARMM) {
      std::cout << "Error: CHARMM-Style parameter is set but EXOTIC-Style parameter file was found.\n"
                "       Either set EXOTIC-Style in config file or change the keyword\n"
                "       \"NONBONDED_MIE\" to \"NONBONDED\" in the parameter files.\n";
    } else {
      std::cout << "Error: EXOTIC-Style parameter is set but CHARMM-Style parameter file was found.\n"
                "       Either set CHARMM-Style in config file or change the keyword\n"
                "       \"NONBONDED\" to \"NONBONDED_MIE\" in the parameter files.\n";
    }
    exit(EXIT_FAILURE);
  }

  //Adjust dih names so lookup finds kind indices rather than term counts
  dih.clean_names();
}

namespace ff_setup
{

std::string FFBase::ReadKind(Reader & param,
                             std::string const& firstKindName)
{
  std::string merged = firstKindName, tmp;
  for (uint k = 1; k < numTerms; k++) {
    param.file >> tmp;
    merged += tmp;
  }
  if (!multi ||
      std::find(name.begin(), name.end(), merged) == name.end())
    name.push_back(merged);
  return merged;
}

std::string FFBase::LoadLine(Reader & param, std::string const& firstVar)
{
  std::string line;
  ReadKind(param, firstVar);
  std::getline(param.file, line);
  return line;
}

void Particle::Read(Reader & param, std::string const& firstVar)
{
  real e, s, e_1_4, s_1_4, dummy1, dummy2;
  uint expN, expN_1_4;
  std::stringstream values(LoadLine(param, firstVar));
  if (isCHARMM()) { //if lj
    values >> dummy1;
  }
  values >> e >> s;
  if (isCHARMM()) {
    expN = ff::part::lj_n;
  } else {
    values >> expN;
  }

  if (isCHARMM()) { //if lj
    values >> dummy2;
  }
  //If undefined in CHARMM, assign 1-4 to full value.
  values >> e_1_4 >> s_1_4;
  if (values.fail()) {
    e_1_4 = e;
    s_1_4 = s;
  }
  values >> expN_1_4;
  if (isCHARMM() || values.fail()) {
    expN_1_4 = expN;
  }
  Add(e, s, expN, e_1_4, s_1_4, expN_1_4);
}

void Particle::Add(real e, real s, const uint expN,
                   real e_1_4, real s_1_4, const uint expN_1_4)
{
  if (isCHARMM()) {
    e *= -1.0;
    s *= RIJ_OVER_2_TO_SIG;
    e_1_4 *= -1.0;
    s_1_4 *= RIJ_OVER_2_TO_SIG;
  }
  epsilon.push_back(EnConvIfCHARMM(e));
  sigma.push_back(s);
  n.push_back(expN);
  epsilon_1_4.push_back(EnConvIfCHARMM(e_1_4));
  sigma_1_4.push_back(s_1_4);
  n_1_4.push_back(expN_1_4);
}

void NBfix::Read(Reader & param, std::string const& firstVar)
{
  real e, s, e_1_4, s_1_4;
#ifdef MIE_INT_ONLY
  uint expN, expN_1_4;
#else
  real expN, expN_1_4;
#endif

  std::stringstream values(LoadLine(param, firstVar));
  values >> e >> s;
  if (isCHARMM()) {
    expN = ff::part::lj_n;
  } else {
    values >> expN;
  }

  values >> e_1_4 >> s_1_4;
  if (values.fail()) {
    e_1_4 = e;
    s_1_4 = s;
  }
  values >> expN_1_4;
  if (isCHARMM() || values.fail()) {
    expN_1_4 = expN;
  }
  Add(e, s, expN, e_1_4, s_1_4, expN_1_4);
}

void NBfix::Add(real e, real s,
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
               )
{
  if (isCHARMM()) {
    e *= -1.0;
    s *= RIJ_TO_SIG;
    e_1_4 *= -1.0;
    s_1_4 *= RIJ_TO_SIG;
  }
  epsilon.push_back(EnConvIfCHARMM(e));
  sigma.push_back(s);
  n.push_back(expN);
  epsilon_1_4.push_back(EnConvIfCHARMM(e_1_4));
  sigma_1_4.push_back(s_1_4);
  n_1_4.push_back(expN_1_4);

}

void Bond::Read(Reader & param, std::string const& firstVar)
{
  real coeff, def;
  ReadKind(param, firstVar);
  param.file >> coeff >> def;
  Add(coeff, def);
}
void Bond::Add(const real coeff, const real def)
{
  fixed.push_back(coeff > FIXED);
  Kb.push_back(EnConvIfCHARMM(coeff));
  b0.push_back(def);
}

void Angle::Read(Reader & param, std::string const& firstVar)
{
  real coeff, def, coeffUB, defUB;
  bool hsUB;
  std::stringstream values(LoadLine(param, firstVar));
  values >> coeff >> def;
  values >> coeffUB >> defUB;

  hsUB = !values.fail();
  Add(coeff, def, hsUB, coeffUB, defUB);
}
void Angle::Add(const real coeff, const real def, const bool hsUB,
                const real coeffUB, const real defUB)
{
  fixed.push_back(coeff > FIXED);
  Ktheta.push_back(EnConvIfCHARMM(coeff));
  theta0.push_back(geom::DegToRad(def));
  hasUB.push_back(hsUB);
  if (hsUB) {
    Kub.push_back(EnConvIfCHARMM(coeffUB));
    bUB0.push_back(defUB);
  } else {
    Kub.push_back(0.0);
    bUB0.push_back(0.0);
  }
}

void Dihedral::Read(Reader & param, std::string const& firstVar)
{
  real coeff, def;
  uint index;
  std::string merged = ReadKind(param, firstVar);
  param.file >> coeff >> index >> def;
  if(index == 0) {
    //set phase shif for n=0 to 90 degree
    // We will have C0 = Kchi (1 + cos(0 * phi + 90)) = Kchi
    def = 90.00;
  }
  Add(merged, coeff, index, def);
  last = merged;
}
void Dihedral::Add(std::string const& merged,
                   const real coeff, const uint index, const real def)
{
  ++countTerms;
  Kchi[merged].push_back(EnConvIfCHARMM(coeff));
  n[merged].push_back(index);
  delta[merged].push_back(geom::DegToRad(def));
}

void Improper::Read(Reader & param, std::string const& firstVar)
{
  real coeff, def;
  std::string merged = ReadKind(param, firstVar);
  //If new value
  if (validname(merged) == false) {
    param.file >> coeff >> def;
    Add(coeff, def);
  }
}
void Improper::Add(const real coeff, const real def)
{
  Komega.push_back(EnConvIfCHARMM(coeff));
  omega0.push_back(def);
}

#ifndef NDEBUG
void Particle::PrintBrief()
{
  std::cout << "\tSigma\t\tEpsilon\t\tN\n";
  std::cout << "#Read\t" << sigma.size() << '\t' << '\t' << epsilon.size()
            << "\t\t" << n.size() << '\n';
  std::cout << "Avg.\t" <<
            std::accumulate(sigma.begin(), sigma.end(), 0.0) / sigma.size()
            << "\t\t" <<
            std::accumulate(epsilon.begin(), epsilon.end(), 0.0) / epsilon.size()
            << "\t\t" <<
            std::accumulate(n.begin(), n.end(), 0.0) / n.size() << "\n\n";
}

void Bond::PrintBrief()
{
  std::cout << "\tKb\t\tb0\n";
  std::cout << "#Read\t" << Kb.size() << '\t' << '\t' << b0.size() << '\n';
  std::cout << "Avg.\t" <<
            std::accumulate(Kb.begin(), Kb.end(), 0.0) / Kb.size() << "\t" <<
            std::accumulate(b0.begin(), b0.end(), 0.0) / b0.size() << "\n\n";
}

void Angle::PrintBrief()
{
  std::cout << "\tKtheta\t\ttheta0\n";
  std::cout << "#Read\t" << Ktheta.size() << '\t' << '\t'
            << theta0.size() << '\n';
  std::cout << "Avg.\t" <<
            std::accumulate(Ktheta.begin(), Ktheta.end(), 0.0) / Ktheta.size()
            << "\t\t" <<
            std::accumulate(theta0.begin(), theta0.end(), 0.0) / theta0.size()
            << "\n\n";
}

void Dihedral::PrintBrief()
{
  std::cout << name.size() << " dihedrals.\n";
}

void Improper::PrintBrief()
{
  std::cout << "\tKomega\t\tomega0\n";
  std::cout << "#Read\t" << Komega.size() << '\t' << '\t'
            << omega0.size() << '\n';
  std::cout << "Avg.\t" <<
            std::accumulate(Komega.begin(),
                            Komega.end(), 0.0) / Komega.size() << "\t" <<
            std::accumulate(omega0.begin(),
                            omega0.end(), 0.0) / omega0.size() << "\n\n";
}
#endif

} //end namespace ff_setup
