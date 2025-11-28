/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "FFSetup.h"

#include <algorithm>
#include <iostream>

#ifndef NDEBUG
#include <numeric>
#endif

#include "GeomLib.h"

// For Exotic style parameter header error checking
#include <regex>

const double EPSILON = 1.0e-4;
const uint FFSetup::CHARMM_ALIAS_IDX = 0;
const uint FFSetup::EXOTIC_ALIAS_IDX = 1;
const std::string FFSetup::paramFileAlias[] = {"CHARMM-Style parameter file",
                                               "Mie-Style parameter file"};
const double ff_setup::KCAL_PER_MOL_TO_K = 503.21959899;
const double ff_setup::RIJ_OVER_2_TO_SIG = 1.7817974362807;
const double ff_setup::RIJ_TO_SIG = 0.890898718;

const double ff_setup::Bond::FIXED = 99999999;
const double ff_setup::Angle::FIXED = 99999999;

// Map variable names to functions
std::map<std::string, ReadableBaseWithFirst *>
FFSetup::SetReadFunctions(const bool isCHARMM) {
  std::map<std::string, ReadableBaseWithFirst *> funct;
  // From CHARMM style file.
  funct["BOND"] = &bond;
  funct["BONDS"] = &bond;
  funct["ANGLE"] = &angle;
  funct["ANGLES"] = &angle;
  funct["DIHEDRAL"] = &dih;
  funct["DIHEDRALS"] = &dih;
  // Unique to Charmm file
  funct["NONBONDED"] = &mie;
  funct["NBFIX"] = &nbfix;
  // Error checking done to ensure isCharmm is true if these are found

  // Unique to exotic file
  funct["NONBONDED_MIE"] = &mie;
  funct["NBFIX_MIE"] = &nbfix;
  // Error checking done to ensure isCharmm is false if these are found

  // Not supported but shouldn't break GOMC
  funct["IMPROPER"] = &imp;
  funct["IMPROPERS"] = &imp;
  funct["CMAP"] = &cmap;
  funct["HBOND"] = &hbond;

  for (sect_it it = funct.begin(); it != funct.end(); ++it) {
    (dynamic_cast<ff_setup::FFBase *>(it->second))->setIsCHARMM(isCHARMM);
  }
  return funct;
}

void FFSetup::Init(const std::vector<config_setup::FileName> &fileName,
                   const bool isCHARMM) {
  sectKind = SetReadFunctions(isCHARMM);
  std::string currSectName = "", varName = "";
  std::string commentChar = "*!";
  std::string commentStr =
      "REMARK ATOM ATOMS MASS set AEXP REXP HAEX AAEX NBOND "
      "CUTNB END CTONN EPS VSWI NBXM INHI";

  for (int p = 0; p < (int)fileName.size(); p++) {
    std::string name = fileName[p].name;

    std::map<std::string, ReadableBaseWithFirst *>::const_iterator sect,
        currSect;

    Reader param(name,
                 paramFileAlias[isCHARMM ? CHARMM_ALIAS_IDX : EXOTIC_ALIAS_IDX],
                 true, &commentStr, true, &commentChar);
    param.open();
    bool shouldImproperWarn = true;
    bool shouldCMapWarn = true;
    bool shouldHBondWarn = true;
    while (param.Read(varName)) {
      sect = sectKind.find(varName);
      if (sect != sectKind.end()) {
        param.SkipLine(); // Skip rest of line for sect. heading
        currSectName = varName;
        currSect = sect; // Save for later calls.
        std::cout << "Reading " << currSectName << " parameters.\n";
        if (shouldImproperWarn &&
            (currSectName == "IMPROPER" || currSectName == "IMPROPERS")) {
          std::cout << "Warning: GOMC does not support IMPROPER!\n";
          shouldImproperWarn = false;
        } else if (shouldCMapWarn && currSectName == "CMAP") {
          std::cout << "Warning: GOMC does not support CMAP!\n";
          shouldCMapWarn = false;
        } else if (shouldHBondWarn && currSectName == "HBOND") {
          std::cout << "Warning: GOMC does not support HBond!\n";
          shouldHBondWarn = false;
        }
        if (isCHARMM) {
          if (hasEnding(currSectName, "MIE")) {
            std::cout << "Error: CHARMM-Style parameter is set but "
                         "Mie-Style parameter header "
                      << currSectName
                      << " was found.\n"
                         "       Either set Mie-Style in config file or "
                         "change the keyword\n"
                         "       "
                      << currSectName << " to "
                      << currSectName.substr(0, currSectName.size() - 4)
                      << " in the parameter files.\n";
            exit(EXIT_FAILURE);
          }
        } else {
          std::regex nbonded("NONBONDED");
          std::regex nbfix("NBFIX");
          if (std::regex_match(currSectName, nbonded) ||
              std::regex_match(currSectName, nbfix)) {
            std::cout << "Error: Mie-Style parameter is set but "
                         "CHARMM-Style parameter header "
                      << currSectName
                      << " was found.\n"
                         "       Either set CHARMM-Style in config file or "
                         "change the keyword\n"
                         "       "
                      << currSectName << " to " << currSectName.append("_MIE")
                      << " in the parameter files.\n";
            exit(EXIT_FAILURE);
          }
        }
      } else {
        if (currSectName != "CMAP" && currSectName != "HBOND")
          currSect->second->Read(param, varName);
      }
    }

    param.close();
  }

  // check if we read nonbonded parameter
  if (mie.sigma.size() == 0) {
    if (isCHARMM) {
      std::cout << "Error: CHARMM-Style parameter is set but Mie-Style "
                   "parameter file was found.\n"
                   "       Either set Mie-Style in config file or change "
                   "the keyword\n"
                   "       \"NONBONDED_MIE\" to \"NONBONDED\" in the parameter "
                   "files.\n";
    } else {
      std::cout << "Error: Mie-Style parameter is set but CHARMM-Style "
                   "parameter file was found.\n"
                   "       Either set CHARMM-Style in config file or change "
                   "the keyword\n"
                   "       \"NONBONDED\" to \"NONBONDED_MIE\" in the parameter "
                   "files.\n";
    }
    exit(EXIT_FAILURE);
  }

  // Adjust dih names so lookup finds kind indices rather than term counts
  dih.clean_names();
}

bool FFSetup::hasEnding(std::string const &fullString,
                        std::string const &ending) {
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare(fullString.length() - ending.length(),
                                    ending.length(), ending));
  } else {
    return false;
  }
}

namespace ff_setup {

std::string FFBase::ReadKind(Reader &param, std::string const &firstKindName) {
  std::string merged = firstKindName, tmp;
  for (uint k = 1; k < numTerms; k++) {
    param.file >> tmp;
    merged += tmp;
  }
  // Skip duplicates unless we allow multiple entries, such as for dihedrals
  if (!multi && std::find(name.begin(), name.end(), merged) != name.end()) {
    std::cout << "Warning: Ignoring duplicate entry for " << merged << " in "
              << param.getFileName().c_str() << ". Using first entry.\n";
  }
  // Might insert duplicates but they will be ignored during execution
  if (!multi || std::find(name.begin(), name.end(), merged) == name.end())
    name.push_back(merged);
  return merged;
}

std::string FFBase::LoadLine(Reader &param, std::string const &firstVar) {
  std::string line;
  ReadKind(param, firstVar);
  std::getline(param.file, line);
  return line;
}

void Particle::Read(Reader &param, std::string const &firstVar) {
  double e, s, e_1_4, s_1_4, dummy1, dummy2;
  double expN, expN_1_4;
  std::stringstream values(LoadLine(param, firstVar));
  if (isCHARMM()) { // if lj
    values >> dummy1;
  }
  values >> e >> s;
  if (isCHARMM()) {
    expN = ff::part::lj_n;
  } else {
    values >> expN;
  }

  if (values.fail()) {
    std::cout << "Error: Incomplete Nonbonded parameters were found in "
                 "parameter file!\n";
    exit(EXIT_FAILURE);
  }

  if (isCHARMM()) { // if lj
    values >> dummy2;
  }
  // If undefined in CHARMM, assign 1-4 to full value.
  values >> e_1_4 >> s_1_4;
  if (values.fail()) {
    e_1_4 = e;
    s_1_4 = s;
  }
  values >> expN_1_4;
  if (isCHARMM() || values.fail()) {
    expN_1_4 = expN;
  }

  // Reset values and perform error checking for consistency.
  // If epsilon = 0 then the sigma and exponent values don't impact the
  // computation. But need the exponents to be large enough that the arithmetic
  // or geometric mean of any pair > 6.0. See FFParticle::Blend() for underlying
  // math.
  double smallVal = 1.0e-20;
  if (std::fabs(e) < smallVal) {
    e = 0.0;
    expN = 12.0; // Set to default (LJ) exponent.
  }
  if (std::fabs(e_1_4) < smallVal) {
    e_1_4 = 0.0;
    expN_1_4 = 12.0; // Set to default (LJ) exponent.
  }
  if ((expN - 6.0) < smallVal) {
    std::cout << "ERROR: Mie exponent must be > 6!\n";
    exit(EXIT_FAILURE);
  }
  if ((expN_1_4 - 6.0) < smallVal) {
    std::cout << "ERROR: Mie exponent for 1-4 interactions must be > 6!\n";
    exit(EXIT_FAILURE);
  }

  Add(e, s, expN, e_1_4, s_1_4, expN_1_4);
}

void Particle::Add(double e, double s, const double expN, double e_1_4,
                   double s_1_4, const double expN_1_4) {
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

void NBfix::Read(Reader &param, std::string const &firstVar) {
  double e, s, e_1_4, s_1_4;
  double expN, expN_1_4;

  std::stringstream values(LoadLine(param, firstVar));
  values >> e >> s;
  if (isCHARMM()) {
    expN = ff::part::lj_n;
  } else {
    values >> expN;
  }

  if (values.fail()) {
    std::cout
        << "Error: Incomplete NBfix parameters were found in parameter file!\n";
    exit(EXIT_FAILURE);
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

void NBfix::Add(double e, double s, const double expN, double e_1_4,
                double s_1_4, const double expN_1_4) {
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

void Bond::Read(Reader &param, std::string const &firstVar) {
  double coeff, def;
  ReadKind(param, firstVar);
  param.file >> coeff >> def;
  if (!param.file.good()) {
    std::cout
        << "Error: Incomplete Bond parameters were found in parameter file!\n";
    exit(EXIT_FAILURE);
  }
  Add(coeff, def);
}
void Bond::Add(const double coeff, const double def) {
  fixed.push_back(coeff > FIXED);
  Kb.push_back(EnConvIfCHARMM(coeff));
  b0.push_back(def);
}

void Angle::Read(Reader &param, std::string const &firstVar) {
  double coeff, def, coeffUB, defUB;
  bool hsUB;
  std::stringstream values(LoadLine(param, firstVar));
  values >> coeff >> def;
  if (values.fail()) {
    std::cout
        << "Error: Incomplete Angle parameters were found in parameter file!\n";
    exit(EXIT_FAILURE);
  }
  values >> coeffUB >> defUB;

  hsUB = !values.fail();
  Add(coeff, def, hsUB, coeffUB, defUB);
}
void Angle::Add(const double coeff, const double def, const bool hsUB,
                const double coeffUB, const double defUB) {
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

void Dihedral::Read(Reader &param, std::string const &firstVar) {
  double coeff, def;
  uint index;
  std::string merged = ReadKind(param, firstVar);
  param.file >> coeff >> index >> def;
  if (!param.file.good()) {
    std::cout << "Error: Incomplete Dihedral parameters were found in "
                 "parameter file!\n";
    exit(EXIT_FAILURE);
  }
  if (index == 0) {
    // set phase shift for n=0 to 90 degree
    // We will have C0 = Kchi (1 + cos(0 * phi + 90)) = Kchi
    // this avoids double counting the C0 (constant offset) term, which is used
    // in force fields like TraPPE
    def = 90.00;
  }
  Add(param.getFileName(), merged, coeff, index, def);
  last = merged;
}
void Dihedral::Add(const std::string &fileName, const std::string &merged,
                   const double coeff, const uint index, const double def) {
  // Check for (and skip) duplicate periodicities for the same dihedral
  // Generate an error and terminate if the duplicate dihedrals have different
  // parameters
  auto Kchi_it = Kchi[merged].begin();
  auto delta_it = delta[merged].begin();
  for (auto it = n[merged].begin(); it != n[merged].end(); ++it) {
    // Found a duplicate dihedral
    if (*it == index) {
      if (std::fabs(*Kchi_it - EnConvIfCHARMM(coeff)) > EPSILON ||
          std::fabs(*delta_it - geom::DegToRad(def)) > EPSILON) {
        std::cout << "Error: Inconsistent dihedral parameters were found in "
                  << fileName.c_str() << " for dihedral " << merged
                  << " with periodicity " << index << "!\n";
        exit(EXIT_FAILURE);
      } else {
        std::cout << "Warning: Ignoring duplicate periodicity of " << index
                  << " for dihedral " << merged << " in " << fileName.c_str()
                  << ".\n";
        return;
      }
    }
    Kchi_it++;
    delta_it++;
  }
  ++countTerms;
  Kchi[merged].push_back(EnConvIfCHARMM(coeff));
  n[merged].push_back(index);
  delta[merged].push_back(geom::DegToRad(def));
}

void Improper::Read(Reader &param, std::string const &firstVar) {
  double coeff, def;
  uint index;
  std::string merged = ReadKind(param, firstVar);
  // If new value
  if (validname(merged)) {
    param.file >> coeff >> index >> def;
    if (!param.file.good()) {
      std::cout << "Error: Incomplete Improper parameters was found in "
                   "parameter file!\n";
      exit(EXIT_FAILURE);
    }
    Add(coeff, def);
  }
}
void Improper::Add(const double coeff, const double def) {
  Komega.push_back(EnConvIfCHARMM(coeff));
  omega0.push_back(def);
}

// Currently dummy method, exact same as improper
void CMap::Read(Reader &param, std::string const &firstVar) {
  double coeff, def;
  uint index;
  std::string merged = ReadKind(param, firstVar);
  // If new value
  if (validname(merged)) {
    param.file >> coeff >> index >> def;
    if (!param.file.good()) {
      std::cout << "Error: Incomplete Improper parameters was found in "
                   "parameter file!\n";
      exit(EXIT_FAILURE);
    }
    Add(coeff, def);
  }
}
// Currently dummy method, exactly the same as improper
void CMap::Add(const double coeff, const double def) {
  Komega.push_back(EnConvIfCHARMM(coeff));
  omega0.push_back(def);
}

// Currently dummy method, exactly the same as improper
void HBond::Read(Reader &param, std::string const &firstVar) {
  double coeff, def;
  uint index;
  std::string merged = ReadKind(param, firstVar);
  // If new value
  if (validname(merged)) {
    param.file >> coeff >> index >> def;
    if (!param.file.good()) {
      std::cout << "Error: Incomplete Improper parameters was found in "
                   "parameter file!\n";
      exit(EXIT_FAILURE);
    }
    Add(coeff, def);
  }
}

// Currently dummy method, exact same as improper
void HBond::Add(const double coeff, const double def) {
  Komega.push_back(EnConvIfCHARMM(coeff));
  omega0.push_back(def);
}

#ifndef NDEBUG
void Particle::PrintBrief() {
  std::cout << "\tSigma\t\tEpsilon\t\tN\n";
  std::cout << "#Read\t" << sigma.size() << '\t' << '\t' << epsilon.size()
            << "\t\t" << n.size() << '\n';
  std::cout << "Avg.\t"
            << std::accumulate(sigma.begin(), sigma.end(), 0.0) / sigma.size()
            << "\t\t"
            << std::accumulate(epsilon.begin(), epsilon.end(), 0.0) /
                   epsilon.size()
            << "\t\t" << std::accumulate(n.begin(), n.end(), 0.0) / n.size()
            << "\n\n";
}

void Bond::PrintBrief() {
  std::cout << "\tKb\t\tb0\n";
  std::cout << "#Read\t" << Kb.size() << '\t' << '\t' << b0.size() << '\n';
  std::cout << "Avg.\t"
            << std::accumulate(Kb.begin(), Kb.end(), 0.0) / Kb.size() << "\t"
            << std::accumulate(b0.begin(), b0.end(), 0.0) / b0.size() << "\n\n";
}

void Angle::PrintBrief() {
  std::cout << "\tKtheta\t\ttheta0\n";
  std::cout << "#Read\t" << Ktheta.size() << '\t' << '\t' << theta0.size()
            << '\n';
  std::cout << "Avg.\t"
            << std::accumulate(Ktheta.begin(), Ktheta.end(), 0.0) /
                   Ktheta.size()
            << "\t\t"
            << std::accumulate(theta0.begin(), theta0.end(), 0.0) /
                   theta0.size()
            << "\n\n";
}

void Dihedral::PrintBrief() { std::cout << name.size() << " dihedrals.\n"; }

void Improper::PrintBrief() {
  std::cout << "\tKomega\t\tomega0\n";
  std::cout << "#Read\t" << Komega.size() << '\t' << '\t' << omega0.size()
            << '\n';
  std::cout << "Avg.\t"
            << std::accumulate(Komega.begin(), Komega.end(), 0.0) /
                   Komega.size()
            << "\t"
            << std::accumulate(omega0.begin(), omega0.end(), 0.0) /
                   omega0.size()
            << "\n\n";
}
#endif

} // end namespace ff_setup
