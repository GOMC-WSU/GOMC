/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FF_TABULATED_H
#define FF_TABULATED_H

#include "BasicTypes.h" //for uint
#include "FFParticle.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

//////////////////////////////////////////////////////////////////////
////////////////////////// Tabulated Potential Style
///////////////////////////////
//////////////////////////////////////////////////////////////////////
// Tabulated potential calculation using NAMD-style file format:
// File format: distance energy force (one line per distance point)
// Energy and virial are computed via linear interpolation from the table
//
// Units:
//   - Distance: Angstroms (Å)
//   - Energy: K (Kelvin) for non-CHARMM, kcal/mol for CHARMM (converted to K)
//   - Force: K/Å for non-CHARMM, kcal/mol/Å for CHARMM (converted to K/Å)
//
// E(r) = interpolated from table
// Vir(r) = F(r) * r, where F(r) is interpolated from table
// F(r) = -dE/dr, so Vir(r) = -r * dE/dr  CalcVir returns Wij = F/r = -dE/dr
//
// When ParaTypeCHARMM=True, energy and force values are converted from
// kcal/mol and kcal/mol/Å to K and K/Å using KCAL_PER_MOL_TO_K constant.

struct TabulatedData {
  std::vector<double> distance; // Distance values
  std::vector<double> energy;   // Energy values
  std::vector<double> force;    // Force values (F = -dE/dr)
  double rMin;                  // Minimum distance in table
  double rMax;                  // Maximum distance in table
  double dr;                    // Spacing between points (if uniform)
  bool uniformSpacing;          // Whether spacing is uniform
  std::string pairType;         // Pair type identifier (e.g., "OO", "CH", etc.)

  TabulatedData()
      : rMin(0.0), rMax(0.0), dr(0.0), uniformSpacing(false), pairType("") {}
};

// Struct-of-Arrays (SoA) layout for better cache locality and GPU-readiness
// This structure stores ALL pair types' data in contiguous flattened arrays
class TabulatedDataSoA {
public:
  // Optimization: Consolidate metadata into a struct (AoS) for better cache
  // locality during index lookup. Keeping data arrays flattened (SoA) for
  // efficient access.
  struct PairParameters {
    double rMin;
    double rMax;
    double dr;
    double invDr;       // Optimization: Pre-calculated 1.0/dr to avoid division
    size_t tableOffset; // Index in flattened arrays
    size_t tableSize;   // Number of points
    bool uniformSpacing;
    // Padding added automatically by compiler for alignment
  };

private:
  // Flattened arrays containing data for ALL pair types
  std::vector<double> allDistances; // All distance tables concatenated
  std::vector<double> allEnergies;  // All energy tables concatenated
  std::vector<double> allForces;    // All force tables concatenated

  // Per-pair consolidated metadata
  std::vector<PairParameters> pairParams;

  // Per-pair type string (kept separate as it's not used in hot loop)
  std::vector<std::string> pairType;

  uint numPairTypes; // Total number of pair types (count²)

public:
  // Constructor
  TabulatedDataSoA(uint count) : numPairTypes(count * count) {
    pairParams.resize(numPairTypes);
    // Initialize with default values
    for (size_t i = 0; i < numPairTypes; ++i) {
      pairParams[i].rMin = 0.0;
      pairParams[i].rMax = 0.0;
      pairParams[i].dr = 0.0;
      pairParams[i].invDr = 0.0;
      pairParams[i].tableOffset = 0;
      pairParams[i].tableSize = 0;
      pairParams[i].uniformSpacing = false;
    }
    pairType.resize(numPairTypes, "");
  }

  // Accessors for backward compatibility
  inline const double *getDistances(uint pairIdx) const {
    return allDistances.data() + pairParams[pairIdx].tableOffset;
  }

  inline const double *getEnergies(uint pairIdx) const {
    return allEnergies.data() + pairParams[pairIdx].tableOffset;
  }

  inline const double *getForces(uint pairIdx) const {
    return allForces.data() + pairParams[pairIdx].tableOffset;
  }

  inline size_t getTableSize(uint pairIdx) const {
    return pairParams[pairIdx].tableSize;
  }

  // New optimized accessor for parameters
  inline const PairParameters &getParams(uint pairIdx) const {
    return pairParams[pairIdx];
  }

  inline double getRMin(uint pairIdx) const { return pairParams[pairIdx].rMin; }
  inline double getRMax(uint pairIdx) const { return pairParams[pairIdx].rMax; }
  inline double getDr(uint pairIdx) const { return pairParams[pairIdx].dr; }
  inline double getInvDr(uint pairIdx) const {
    return pairParams[pairIdx].invDr;
  }
  inline bool isUniformSpacing(uint pairIdx) const {
    return pairParams[pairIdx].uniformSpacing;
  }
  inline const std::string &getPairType(uint pairIdx) const {
    return pairType[pairIdx];
  }

  // Setter for initialization
  void setPairData(uint pairIdx, const std::vector<double> &dist,
                   const std::vector<double> &energy,
                   const std::vector<double> &force, double rMinVal,
                   double rMaxVal, double drVal, bool uniformSpace,
                   const std::string &pairTypeVal) {
    // Set metadata
    pairParams[pairIdx].rMin = rMinVal;
    pairParams[pairIdx].rMax = rMaxVal;
    pairParams[pairIdx].dr = drVal;
    pairParams[pairIdx].invDr = (drVal > 1e-10) ? (1.0 / drVal) : 0.0;
    pairParams[pairIdx].uniformSpacing = uniformSpace;
    pairType[pairIdx] = pairTypeVal;

    // Set offset and size
    pairParams[pairIdx].tableOffset = allDistances.size();
    pairParams[pairIdx].tableSize = dist.size();

    // Append data to flattened arrays
    allDistances.insert(allDistances.end(), dist.begin(), dist.end());
    allEnergies.insert(allEnergies.end(), energy.begin(), energy.end());
    allForces.insert(allForces.end(), force.begin(), force.end());
  }
};

struct FF_TABULATED : public FFParticle {
public:
  FF_TABULATED(Forcefield &ff)
      : FFParticle(ff), tableData(NULL), tableData_1_4(NULL),
        nbtableData(NULL) {}

  virtual ~FF_TABULATED() {
    delete tableData;
    delete tableData_1_4;
  }

  virtual void Init(ff_setup::Particle const &mie,
                    ff_setup::NBfix const &nbfix);

  // Initialize with NBtable data available for lookup
  void InitWithNBtable(ff_setup::Particle const &mie,
                       ff_setup::NBfix const &nbfix,
                       const ff_setup::NBtable &nbtable) {
    // Set the NBtable data first, before calling Init
    nbtableData = &nbtable;
    // Now call the regular Init
    Init(mie, nbfix);
  }

  // Build mapping from (kind1, kind2) pairs to pair type names from NBtable
  // Uses the provided Mie particle name list to resolve atom type strings to
  // kind indices. This avoids assuming NBtable entries are in the same
  // sequential order as the internal FlatIndex ordering.
  // OPTIMIZATION: Uses hash map for O(1) lookup instead of O(n) nested loop
  void BuildPairTypeMap(const ff_setup::Particle &mie) {
    pairTypeMap.clear();

    if (nbtableData == NULL) {
      throw std::runtime_error(
          "Fatal Error: FF_TABULATED initialized without NBtable data. "
          "Ensure InitWithNBtable is called instead of Init.");
    }

    std::cout << "Building Pair Type Map for Tabulated Potentials:"
              << std::endl;
    std::cout << "Total particle kinds: " << count << std::endl;
    std::cout << "Total NBtable entries: " << nbtableData->GetPairCount()
              << std::endl;

    // Build reverse lookup map for O(1) access - OPTIMIZATION
    std::map<std::string, uint> nameToKind;
    for (uint k = 0; k < mie.getnamecnt(); ++k) {
      nameToKind[mie.getname(k)] = k;
      std::cout << "'" << mie.getname(k) << "'";
      if (k + 1 < mie.getnamecnt())
        std::cout << ", ";
    }
    std::cout << std::endl;

    if (nbtableData->GetPairCount() == 0) {
      throw std::runtime_error("Fatal Error: NBTable has no entries. Check "
                               "the tabulated potential file.");
    }

    // Now iterate with O(1) lookups instead of O(n²)
    for (size_t t = 0; t < nbtableData->GetPairCount(); ++t) {
      std::string atom1 = nbtableData->GetAtomType1(t);
      std::string atom2 = nbtableData->GetAtomType2(t);
      std::string pairType = nbtableData->GetTablePairName(t);

      std::cout << "NBtable entry " << t << ": atoms=('" << atom1 << "','"
                << atom2 << "') -> pairType='" << pairType << "'" << std::endl;

      // Use hash map lookup instead of nested loop - O(1) instead of O(n)
      std::map<std::string, uint>::const_iterator it1 = nameToKind.find(atom1);
      std::map<std::string, uint>::const_iterator it2 = nameToKind.find(atom2);

      if (it1 == nameToKind.end() || it2 == nameToKind.end()) {
        std::cout << "  (Available mie names: ";
        for (uint k = 0; k < mie.getnamecnt(); ++k) {
          std::cout << "'" << mie.getname(k) << "'";
          if (k + 1 < mie.getnamecnt())
            std::cout << ", ";
        }
        std::cout << ")" << std::endl;
        throw std::runtime_error(
            "Fatal Error: Could not find atom types '" + atom1 + "' and '" +
            atom2 + "' in particle type list (pairType: '" + pairType + "')");
      }

      uint kindA = it1->second;
      uint kindB = it2->second;
      uint minKind = std::min(kindA, kindB);
      uint maxKind = std::max(kindA, kindB);

      std::stringstream key;
      key << minKind << "_" << maxKind;
      std::string keyStr = key.str();

      pairTypeMap[keyStr] = pairType;
    }

    std::cout << "Map Building Complete: " << pairTypeMap.size()
              << " entries in map" << std::endl;
    std::cout << std::endl;
  }

  // Method to register tabulated potential file for a specific pair type
  void RegisterPairTableFile(const uint kind1, const uint kind2,
                             const std::string &filename) {
    uint minKind = std::min(kind1, kind2);
    uint maxKind = std::max(kind1, kind2);

    std::stringstream ss;
    ss << minKind << "_" << maxKind;
    pairTableFiles[ss.str()] = filename;
  }

  // Method to register tabulated potential files from NBtable setup data
  void RegisterPairTableFiles(const ff_setup::NBtable &nbtable) {
    // Extract pair information from nbtable
    // The nbtable contains a list of pair names in format "type1_type2"
    for (size_t i = 0; i < nbtable.getnamecnt(); i++) {
      std::string pairName = nbtable.getname(i);
      // Parse "type1_type2" format
      size_t underscore = pairName.find('_');
      if (underscore != std::string::npos) {
        std::string kind1Str = pairName.substr(0, underscore);
        std::string kind2Str = pairName.substr(underscore + 1);

        try {
          uint kind1 = std::stoul(kind1Str);
          uint kind2 = std::stoul(kind2Str);

          std::cout << "Registered pair type (" << kind1 << ", " << kind2 << ")"
                    << std::endl;
        } catch (const std::exception &e) {
          std::cout << "Warning: Could not parse pair type from NBtable entry: "
                    << pairName << std::endl;
        }
      }
    }
  }

  virtual double CalcEn(const double distSq, const uint kind1, const uint kind2,
                        const double lambda) const;
  virtual double CalcVir(const double distSq, const uint kind1,
                         const uint kind2, const double lambda) const;
  virtual void CalcAdd_1_4(double &en, const double distSq, const uint kind1,
                           const uint kind2) const;

  // coulomb interaction functions
  virtual double CalcCoulomb(const double distSq, const uint kind1,
                             const uint kind2, const double qi_qj_Fact,
                             const double lambda, const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const uint kind1,
                                const uint kind2, const double qi_qj,
                                const double lambda, const uint b) const;
  virtual void CalcCoulombAdd_1_4(double &en, const double distSq,
                                  const double qi_qj_Fact, const bool NB) const;

  //! Returns Ezero, no energy correction
  virtual double EnergyLRC(const uint kind1, const uint kind2) const {
    return 0.0;
  }
  //!!Returns Ezero, no virial correction
  virtual double VirialLRC(const uint kind1, const uint kind2) const {
    return 0.0;
  }
  //! Returns zero for impulse pressure correction term for a kind pair
  virtual double ImpulsePressureCorrection(const uint kind1,
                                           const uint kind2) const {
    return 0.0;
  }

  // Calculate the dE/dlambda for vdw energy
  virtual double CalcdEndL(const double distSq, const uint kind1,
                           const uint kind2, const double lambda) const;
  // Calculate the dE/dlambda for Coulomb energy
  virtual double CalcCoulombdEndL(const double distSq, const uint kind1,
                                  const uint kind2, const double qi_qj_Fact,
                                  const double lambda, uint b) const;

protected:
  virtual double CalcEn(const double distSq, const uint index) const;
  virtual double CalcVir(const double distSq, const uint index) const;
  virtual double CalcCoulomb(const double distSq, const double qi_qj_Fact,
                             const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const double qi_qj,
                                uint b) const;

  // Helper functions for tabulated potential
  void ReadTabulatedFile(const std::string &filename, uint pairIndex,
                         const std::string &pairTypeStr, const bool isCHARMM);
  double InterpolateEnergy(uint pairIndex, const double dist) const;
  double InterpolateForce(uint pairIndex, const double dist) const;
  uint FindTableIndex(uint pairIndex, const double dist) const;
  std::string GetPairTypeString(const uint kind1, const uint kind2) const;

  // Storage for tabulated data
  TabulatedDataSoA *tableData;     // Tabulated data for each pair type
  TabulatedDataSoA *tableData_1_4; // Tabulated data for 1-4 interactions

  // Map to store per-pair tabulated potential file names
  // Key format: "type1_type2" where type1 <= type2
  std::map<std::string, std::string> pairTableFiles;

  // Reference to NBtable data for pair type lookups
  const ff_setup::NBtable *nbtableData;

  // Map from "(kind1_kind2)" to pair type name from NBtable
  std::map<std::string, std::string> pairTypeMap;
};

/**
 * @brief Read tabulated potential data from a file
 * @param filename Name of the file to read
 * @param pairIndex Index of the pair in the SoA structure
 * @param pairTypeStr String identifier for the pair type
 * @param isCHARMM Whether the data is in CHARMM format
 * @note This function reads all tabulated potential data from a file (even if
 * not needed).
 * @note The data is stored in the SoA structure via setPairData.
 */
inline void FF_TABULATED::ReadTabulatedFile(const std::string &filename,
                                            uint pairIndex,
                                            const std::string &pairTypeStr,
                                            const bool isCHARMM) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    std::cout << "ERROR: Cannot open tabulated potential file: " << filename
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  double r, e, f;
  bool foundPair = false;
  int dbgPrinted = 0;

  // Temporary vectors to hold data before storing in SoA
  std::vector<double> tempDistance;
  std::vector<double> tempEnergy;
  std::vector<double> tempForce;

  // Unit conversion factor for CHARMM format
  // Energy: kcal/mol -> K, Force: kcal/mol/Å -> K/Å
  double energyConversion = 1.0;
  double forceConversion = 1.0;
  if (isCHARMM) {
    energyConversion = ff_setup::KCAL_PER_MOL_TO_K;
    forceConversion = ff_setup::KCAL_PER_MOL_TO_K;
    std::cout
        << "DEBUG: Using CHARMM unit conversions for tabulated potential file: "
        << filename << std::endl;
  }

  // Read file line by line looking for the specified pair type
  while (std::getline(file, line)) {
    // Trim leading whitespace
    size_t start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
      continue; // Empty line

    line = line.substr(start);

    // Check for TYPE keyword
    if (line.substr(0, 4) == "TYPE") {
      // If we were already reading data and now found a new TYPE, stop
      if (foundPair && !tempDistance.empty()) {
        break;
      }

      // Parse the pair type from the line
      std::istringstream typeStream(line);
      std::string typeKeyword, typeStr;
      typeStream >> typeKeyword >> typeStr;

      // std::cout << "DEBUG: Found TYPE keyword with pair: '" << typeStr << "'"
      //           << " (looking for '" << pairTypeStr << "')" << std::endl;

      // Check if this is the pair type we're looking for
      if (typeStr == pairTypeStr) {
        foundPair = true;
        //  std::cout << "Found pair type: " << typeStr << std::endl;
        continue;
      } else if (foundPair) {
        // We were reading this pair and now found a different one, so stop
        break;
      }
      continue;
    }

    // Skip comments (lines starting with # or !)
    if (line[0] == '#' || line[0] == '!')
      continue;

    // If we haven't found our pair type yet, skip this line
    if (!foundPair)
      continue;

    // Try to read three numbers: distance, energy, force
    std::istringstream iss(line);
    if (iss >> r >> e >> f) {
      tempDistance.push_back(r);
      // Convert energy from kcal/mol to K if CHARMM format
      tempEnergy.push_back(e * energyConversion);
      // Convert force from kcal/mol/Å to K/Å if CHARMM format
      tempForce.push_back(f * forceConversion);
    }
  }

  file.close();

  if (!foundPair) {
    std::cout << "ERROR: Pair type '" << pairTypeStr << "' not found in "
              << "tabulated potential file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  if (tempDistance.empty()) {
    std::cout << "ERROR: No valid data found for pair type '" << pairTypeStr
              << "' in tabulated potential file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  // Find min and max distances using OpenMP parallelization
  // OPTIMIZATION: Parallelize min/max finding for large tables
  double localMin = std::numeric_limits<double>::max();
  double localMax = std::numeric_limits<double>::lowest();

#ifdef _OPENMP
#pragma omp parallel for reduction(min : localMin) reduction(max : localMax)
#endif
  for (size_t i = 0; i < tempDistance.size(); i++) {
    if (tempDistance[i] < localMin)
      localMin = tempDistance[i];
    if (tempDistance[i] > localMax)
      localMax = tempDistance[i];
  }

  double rMin = localMin;
  double rMax = localMax;

  // Check if spacing is uniform
  // OPTIMIZATION: Use SIMD-friendly loop for uniform spacing check
  double dr = 0.0;
  bool uniformSpacing = false;
  if (tempDistance.size() > 1) {
    dr = tempDistance[1] - tempDistance[0];
    uniformSpacing = true;

    // Note: Attempted SIMD optimization here, but conditional logic prevents
    // vectorization
    for (size_t i = 2; i < tempDistance.size(); i++) {
      double expected = tempDistance[0] + i * dr;
      if (std::fabs(tempDistance[i] - expected) > 1.0e-6) {
        uniformSpacing = false;
        break; // Early exit when non-uniform spacing detected
      }
    }
  }

  std::string unitStr = isCHARMM ? "kcal/mol (converted to K)" : "K";
  std::cout << "Info: Loaded tabulated potential for pair type '" << pairTypeStr
            << "' from " << filename << " with " << tempDistance.size()
            << " points, range [" << rMin << ", " << rMax
            << "] Angstrom, energy units: " << unitStr << std::endl;

  // Store data in SoA structure
  tableData->setPairData(pairIndex, tempDistance, tempEnergy, tempForce, rMin,
                         rMax, dr, uniformSpacing, pairTypeStr);
}

inline uint FF_TABULATED::FindTableIndex(uint pairIndex,
                                         const double dist) const {
  // Optimization: Access all metadata from a single struct to minimize cache
  // misses
  const TabulatedDataSoA::PairParameters &params =
      tableData->getParams(pairIndex);

  if (dist <= params.rMin)
    return 0;
  if (dist >= params.rMax)
    return params.tableSize - 2; // Return second-to-last for interpolation

  // Binary search for non-uniform spacing, or direct calculation for uniform
  if (params.uniformSpacing) {
    // Optimization: Use multiplication by inverse dr instead of division
    uint idx = static_cast<uint>((dist - params.rMin) * params.invDr);
    return std::min(idx, static_cast<uint>(params.tableSize - 2));
  } else {
    // Binary search
    const double *distances = tableData->getDistances(pairIndex);
    uint left = 0;
    uint right = params.tableSize - 1;
    while (right - left > 1) {
      uint mid = (left + right) / 2;
      if (distances[mid] <= dist)
        left = mid;
      else
        right = mid;
    }
    return left;
  }
}

inline double FF_TABULATED::InterpolateEnergy(uint pairIndex,
                                              const double dist) const {
  const TabulatedDataSoA::PairParameters &params =
      tableData->getParams(pairIndex);

  // Safety check: if table is empty, return zero energy
  if (params.tableSize == 0) {
    throw std::runtime_error(
        "Fatal Error: InterpolateEnergy called with empty "
        "table for pair type '" +
        tableData->getPairType(pairIndex) +
        "'. Check force field file for missing NBTable entry.");
  }

  const double *distances = tableData->getDistances(pairIndex);
  const double *energies = tableData->getEnergies(pairIndex);

  if (dist <= params.rMin)
    return energies[0];
  if (dist >= params.rMax)
    return energies[params.tableSize - 1];

  uint idx = FindTableIndex(pairIndex, dist);
  double r1 = distances[idx];
  double r2 = distances[idx + 1];
  double e1 = energies[idx];
  double e2 = energies[idx + 1];

  // Check interpolation type
  if (forcefield.interpolationType == "cubic" && idx >= 1 &&
      idx + 2 < params.tableSize) {
    // Cubic interpolation using 4 points: idx-1, idx, idx+1, idx+2
    double x0 = distances[idx - 1];
    double x1 = r1;
    double x2 = r2;
    double x3 = distances[idx + 2];
    double y0 = energies[idx - 1];
    double y1 = e1;
    double y2 = e2;
    double y3 = energies[idx + 2];

    // Lagrange interpolation
    double term0 = y0 * (dist - x1) * (dist - x2) * (dist - x3) /
                   ((x0 - x1) * (x0 - x2) * (x0 - x3));
    double term1 = y1 * (dist - x0) * (dist - x2) * (dist - x3) /
                   ((x1 - x0) * (x1 - x2) * (x1 - x3));
    double term2 = y2 * (dist - x0) * (dist - x1) * (dist - x3) /
                   ((x2 - x0) * (x2 - x1) * (x2 - x3));
    double term3 = y3 * (dist - x0) * (dist - x1) * (dist - x2) /
                   ((x3 - x0) * (x3 - x1) * (x3 - x2));

    return term0 + term1 + term2 + term3;
  } else {
    // Linear interpolation
    double t = (dist - r1) / (r2 - r1);
    return e1 + t * (e2 - e1);
  }
}

inline double FF_TABULATED::InterpolateForce(uint pairIndex,
                                             const double dist) const {
  const TabulatedDataSoA::PairParameters &params =
      tableData->getParams(pairIndex);

  // Safety check: if table is empty, return zero force
  if (params.tableSize == 0) {
    std::cout
        << "WARNING: InterpolateForce called with empty table for pair type '"
        << tableData->getPairType(pairIndex) << "'" << std::endl;
    return 0.0;
  }

  const double *distances = tableData->getDistances(pairIndex);
  const double *forces = tableData->getForces(pairIndex);

  if (dist <= params.rMin)
    return forces[0];
  if (dist >= params.rMax)
    return forces[params.tableSize - 1];

  uint idx = FindTableIndex(pairIndex, dist);
  double r1 = distances[idx];
  double r2 = distances[idx + 1];
  double f1 = forces[idx];
  double f2 = forces[idx + 1];

  // Check interpolation type
  if (forcefield.interpolationType == "cubic" && idx >= 1 &&
      idx + 2 < params.tableSize) {
    // Cubic interpolation using 4 points: idx-1, idx, idx+1, idx+2
    double x0 = distances[idx - 1];
    double x1 = r1;
    double x2 = r2;
    double x3 = distances[idx + 2];
    double y0 = forces[idx - 1];
    double y1 = f1;
    double y2 = f2;
    double y3 = forces[idx + 2];

    // Lagrange interpolation
    double term0 = y0 * (dist - x1) * (dist - x2) * (dist - x3) /
                   ((x0 - x1) * (x0 - x2) * (x0 - x3));
    double term1 = y1 * (dist - x0) * (dist - x2) * (dist - x3) /
                   ((x1 - x0) * (x1 - x2) * (x1 - x3));
    double term2 = y2 * (dist - x0) * (dist - x1) * (dist - x3) /
                   ((x2 - x0) * (x2 - x1) * (x2 - x3));
    double term3 = y3 * (dist - x0) * (dist - x1) * (dist - x2) /
                   ((x3 - x0) * (x3 - x1) * (x3 - x2));

    return term0 + term1 + term2 + term3;
  } else {
    // Linear interpolation
    double t = (dist - r1) / (r2 - r1);
    return f1 + t * (f2 - f1);
  }
}

inline void FF_TABULATED::Init(ff_setup::Particle const &mie,
                               ff_setup::NBfix const &nbfix) {
  // Initialize sigma and epsilon (needed for compatibility)
  FFParticle::Init(mie, nbfix);

  // Build the pair type map NOW that count is set and nbtableData is available
  if (nbtableData != NULL) {
    BuildPairTypeMap(mie);
  } else {
    std::cout << "WARNING: nbtableData is NULL in Init()" << std::endl;
  }

  // NEW: Allocate SoA structures
  tableData = new TabulatedDataSoA(count);
  tableData_1_4 = new TabulatedDataSoA(count);

  // Check if tabulated potential file is specified
  if (forcefield.tabulatedPotentialFile.empty()) {
    std::cout << "ERROR: Tabulated potential file not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Initialize all tables with data from the appropriate files
  // For each pair type (i, j), get the corresponding table from the
  // consolidated file Note: We iterate ALL pairs (i,j) and use the symmetric
  // mapping since FlatIndex(i,j) can be called with any order during runtime.
  // OPTIMIZATION: Parallelize with OpenMP for multicore performance
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
  for (uint i = 0; i < count; i++) {
    for (uint j = 0; j < count; j++) {
      uint index = FlatIndex(i, j);

      // Check if this pair exists in the pairTypeMap using min/max ordering
      uint minKind = std::min(i, j);
      uint maxKind = std::max(i, j);
      std::stringstream keyStream;
      keyStream << minKind << "_" << maxKind;
      std::string key = keyStream.str();

      std::map<std::string, std::string>::const_iterator it =
          pairTypeMap.find(key);
      if (it == pairTypeMap.end()) {
#ifdef _OPENMP
#pragma omp critical(console_output)
#endif
        {
          std::cout << "Warning: Pair (" << i << ", " << j << ") [key '" << key
                    << "'] not found in pairTypeMap; skipping table load"
                    << std::endl;
        }
        // For unknown pairs, set empty data (will be handled by setPairData
        // with empty vectors)
        std::vector<double> empty;
        std::string unknownPairType =
            "UnknownPair: " + mie.getname(i) + "_" + mie.getname(j);
#ifdef _OPENMP
#pragma omp critical(file_reading)
#endif
        {
          tableData->setPairData(index, empty, empty, empty, 0.0, 0.0, 0.0,
                                 false, unknownPairType);
          tableData_1_4->setPairData(index, empty, empty, empty, 0.0, 0.0, 0.0,
                                     false, unknownPairType);
        }
        continue;
      }

      // Construct the pair type string (e.g., "OO", "CH", etc.)
      // This should match the TYPE keyword in the tabulated file
      std::string pairType = it->second;

      // Read the tabulated potential data for this pair from the consolidated
      // file Only read once (check if already loaded to avoid re-reading)
      // Critical section to ensure thread-safe file I/O
#ifdef _OPENMP
#pragma omp critical(file_reading)
#endif
      {
        if (tableData->getTableSize(index) == 0) {
          ReadTabulatedFile(forcefield.tabulatedPotentialFile, index, pairType,
                            forcefield.isCHARMM);
          // Copy to 1-4 table (same data for now)
          tableData_1_4->setPairData(
              index,
              std::vector<double>(tableData->getDistances(index),
                                  tableData->getDistances(index) +
                                      tableData->getTableSize(index)),
              std::vector<double>(tableData->getEnergies(index),
                                  tableData->getEnergies(index) +
                                      tableData->getTableSize(index)),
              std::vector<double>(tableData->getForces(index),
                                  tableData->getForces(index) +
                                      tableData->getTableSize(index)),
              tableData->getRMin(index), tableData->getRMax(index),
              tableData->getDr(index), tableData->isUniformSpacing(index),
              pairType);
        }
      }
    }
  }
}

// Helper function to get the pair type string from atom type indices
inline std::string FF_TABULATED::GetPairTypeString(const uint kind1,
                                                   const uint kind2) const {
  // Normalize the pair (use min, max order)
  uint minKind = std::min(kind1, kind2);
  uint maxKind = std::max(kind1, kind2);

  std::stringstream keyStream;
  keyStream << minKind << "_" << maxKind;
  std::string key = keyStream.str();

  std::cout << "DEBUG GetPairTypeString: Looking for key '" << key
            << "' in map with " << pairTypeMap.size() << " entries"
            << std::endl;

  // Look up in the pre-built map
  std::map<std::string, std::string>::const_iterator it = pairTypeMap.find(key);
  if (it != pairTypeMap.end()) {
    std::cout << "  FOUND: Pair (" << kind1 << ", " << kind2 << ") -> '"
              << it->second << "'" << std::endl;
    return it->second;
  }

  // Fallback: return indices as string if map lookup fails
  std::cout << "  NOT FOUND: Pair (" << kind1 << ", " << kind2
            << ") not in map, using fallback" << std::endl;

  // Print map contents for debugging
  if (pairTypeMap.empty()) {
    std::cout << "  Map is EMPTY!" << std::endl;
  } else {
    std::cout << "  Map contains:" << std::endl;
    for (auto &entry : pairTypeMap) {
      std::cout << "    '" << entry.first << "' -> '" << entry.second << "'"
                << std::endl;
    }
  }

  return key;
}

inline double FF_TABULATED::CalcEn(const double distSq, const uint kind1,
                                   const uint kind2,
                                   const double lambda) const {
  if (forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  if (lambda >= 0.999999) {
    return CalcEn(distSq, index);
  }

  // For lambda scaling, use soft-core approach
  // similar to other force fields
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef =
      forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double dist = sqrt(softRsq);

  double en = InterpolateEnergy(index, dist);
  return lambda * en;
}

inline double FF_TABULATED::CalcEn(const double distSq,
                                   const uint index) const {
  double dist = sqrt(distSq);
  return InterpolateEnergy(index, dist);
}

inline double FF_TABULATED::CalcVir(const double distSq, const uint kind1,
                                    const uint kind2,
                                    const double lambda) const {
  if (forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  if (lambda >= 0.999999) {
    return CalcVir(distSq, index);
  }

  // For lambda scaling
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef =
      forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double dist = sqrt(softRsq);
  double correction = distSq / softRsq;

  double force = InterpolateForce(index, dist);
  // Use the same soft-core scaling as other
  // forcefields: scale the per-index virial (which
  // is F/r) by correction^2 and lambda
  double vir = lambda * correction * correction * CalcVir(softRsq, index);
  return vir;
}

inline double FF_TABULATED::CalcVir(const double distSq,
                                    const uint index) const {
  double dist = sqrt(distSq);
  double force = InterpolateForce(index, dist);
  // CalcVir should return the scalar Wij =
  // (-dE/dr)/r so that r_vector * Wij = force
  // vector. Since InterpolateForce returns F =
  // -dE/dr, we return Wij = F / r.
  if (dist > 0.0)
    return force / dist;
  else
    return 0.0;
}

inline void FF_TABULATED::CalcAdd_1_4(double &en, const double distSq,
                                      const uint kind1,
                                      const uint kind2) const {
  if (forcefield.rCutSq < distSq)
    return;

  uint index = FlatIndex(kind1, kind2);
  double dist = sqrt(distSq);
  en += InterpolateEnergy(index, dist);
}

inline double FF_TABULATED::CalcCoulomb(const double distSq, const uint kind1,
                                        const uint kind2,
                                        const double qi_qj_Fact,
                                        const double lambda,
                                        const uint b) const {
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if (lambda >= 0.999999) {
    return CalcCoulomb(distSq, qi_qj_Fact, b);
  }

  double en = 0.0;
  if (forcefield.sc_coul) {
    uint index = FlatIndex(kind1, kind2);
    double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
    sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    en = lambda * CalcCoulomb(softRsq, qi_qj_Fact, b);
  } else {
    en = lambda * CalcCoulomb(distSq, qi_qj_Fact, b);
  }
  return en;
}

inline double FF_TABULATED::CalcCoulomb(const double distSq,
                                        const double qi_qj_Fact,
                                        const uint b) const {
  if (forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return qi_qj_Fact * erfc(val) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_Fact / dist;
  }
}

inline double FF_TABULATED::CalcCoulombVir(const double distSq,
                                           const uint kind1, const uint kind2,
                                           const double qi_qj,
                                           const double lambda,
                                           const uint b) const {
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if (lambda >= 0.999999) {
    return CalcCoulombVir(distSq, qi_qj, b);
  }

  double vir = 0.0;
  if (forcefield.sc_coul) {
    uint index = FlatIndex(kind1, kind2);
    double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
    sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    vir = lambda * correction * correction * CalcCoulombVir(softRsq, qi_qj, b);
  } else {
    vir = lambda * CalcCoulombVir(distSq, qi_qj, b);
  }
  return vir;
}

inline double FF_TABULATED::CalcCoulombVir(const double distSq,
                                           const double qi_qj,
                                           const uint b) const {
  if (forcefield.ewald) {
    double dist = sqrt(distSq);
    double constValue = forcefield.alpha[b] * M_2_SQRTPI;
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = erfc(forcefield.alpha[b] * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    double dist = sqrt(distSq);
    return qi_qj / (distSq * dist);
  }
}

inline void FF_TABULATED::CalcCoulombAdd_1_4(double &en, const double distSq,
                                             const double qi_qj_Fact,
                                             const bool NB) const {
  if (forcefield.rCutSq < distSq)
    return;

  double dist = sqrt(distSq);
  if (NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}

inline double FF_TABULATED::CalcdEndL(const double distSq, const uint kind1,
                                      const uint kind2,
                                      const double lambda) const {
  if (forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef =
      forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double dist = sqrt(softRsq);

  double fCoef = lambda * forcefield.sc_alpha * forcefield.sc_power / 6.0;
  fCoef *= pow(1.0 - lambda, forcefield.sc_power - 1.0) * sigma6 /
           (softRsq * softRsq);

  double en = InterpolateEnergy(index, dist);
  double force = InterpolateForce(index, dist);
  // CalcVir(index) returns Wij = F/r, so use vir =
  // Wij
  double vir = 0.0;
  if (dist > 0.0)
    vir = force / dist;

  double dhdl = en + fCoef * vir;
  return dhdl;
}

inline double FF_TABULATED::CalcCoulombdEndL(const double distSq,
                                             const uint kind1, const uint kind2,
                                             const double qi_qj_Fact,
                                             const double lambda,
                                             uint b) const {
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  double dhdl = 0.0;
  if (forcefield.sc_coul) {
    uint index = FlatIndex(kind1, kind2);
    double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
    sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double fCoef = lambda * forcefield.sc_alpha * forcefield.sc_power / 6.0;
    fCoef *= pow(1.0 - lambda, forcefield.sc_power - 1.0) * sigma6 /
             (softRsq * softRsq);
    dhdl = CalcCoulomb(softRsq, qi_qj_Fact, b) +
           fCoef * CalcCoulombVir(softRsq, qi_qj_Fact, b);
  } else {
    dhdl = CalcCoulomb(distSq, qi_qj_Fact, b);
  }
  return dhdl;
}

#endif /*FF_TABULATED_H*/
