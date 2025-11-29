/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef PRNG_H
#define PRNG_H

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "BasicTypes.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "EnsemblePreprocessor.h"
#include "MersenneTwister.h"
#include "MoleculeLookup.h"
#include "MoveConst.h"
#include "XYZArray.h"

// Wrapper class for our random numbers
//
// NOTE: gcc is not generating correct code when we call a Mersenne Twister
// function inside a function call. In other words, when the return value is
// used directly as a function parameter. So, many of the functions were
// rewritten to assign the result to a variable. This looks cumbersome, but it
// at least produces correct results. This is the case with some older versions
// of gcc, such as gcc 7.3. It's possible that newer versions of gcc don't have
// this issue, but we need to support older compilers, too, so...
//
class PRNG {
public:
  PRNG(MoleculeLookup &molLook)
      : molLookRef(molLook), hasSecondGaussian(false) {}
  ~PRNG(void) { delete gen; }

  void Init(MTRand *prng) { gen = prng; }

  // Saves the current state of the PRNG as ./filename
  void saveState(const char *filename);

  ////
  // BASIC GENERATION
  ////

  // Standard double generation on [0,1.0]
  double operator()() { return (*gen)(); }

  // Generate a double on [0,bound]
  double rand(const double bound) { return gen->rand(bound); }

  // Generate a double on [0,bound)
  double randExc(const double bound) { return gen->randExc(bound); }

  // Generate an unsigned int on [0,bound]
  uint randInt(const uint bound) { return (uint)(gen->randInt(bound)); }

  // Generate an unsigned int on [0,bound)
  uint randIntExc(const uint bound) { return (uint)(gen->randInt(bound - 1)); }

  // Generates a double on [-bound,bound]
  double Sym(const double bound) { return gen->rand(2.0 * bound) - bound; }

  // return between (-bound, bound)
  double SymExc(double bound) { return gen->randExc(2.0 * bound) - bound; }

  /////////////////////////////
  //   GENERATION FUNCTIONS  //
  /////////////////////////////

  XYZ SymXYZ(const double bound) {
    double bound2 = 2.0 * bound;
    double x = gen->rand(bound2) - bound;
    double y = gen->rand(bound2) - bound;
    double z = gen->rand(bound2) - bound;
    return XYZ(x, y, z);
  }

  // Used to pick first position of
  void FillWithRandom(XYZArray &loc, const uint len, BoxDimensions const &dims,
                      const uint b) {
    for (uint i = 0; i < len; ++i) {
      double x = randExc(dims.axis.x[b]);
      double y = randExc(dims.axis.y[b]);
      double z = randExc(dims.axis.z[b]);
      loc.Set(i, x, y, z);
      loc.Set(i, dims.TransformSlant(loc.Get(i), b));
    }
  }

  void FillWithRandom(XYZ &loc, BoxDimensions const &dims, const uint b) {
    double x = randExc(dims.axis.x[b]);
    double y = randExc(dims.axis.y[b]);
    double z = randExc(dims.axis.z[b]);
    XYZ temp(x, y, z);
    loc = dims.TransformSlant(temp, b);
  }

  void FillWithRandomInCavity(XYZ &loc, XYZ const &cavDim) {
    double x = Sym(cavDim.x / 2.0);
    double y = Sym(cavDim.y / 2.0);
    double z = Sym(cavDim.z / 2.0);
    XYZ temp(x, y, z);
    loc = temp;
  }

  // generate random coordinate in cavity with dimension of cavDim and center
  // of cavCenter
  void FillWithRandomInCavity(XYZArray &loc, const uint len, XYZ const &cavDim,
                              XYZ const &cavCenter) {
    // generate random trial in range of [-cavDim/2, +cavDim/2]
    for (uint i = 0; i < len; ++i) {
      double x = Sym(cavDim.x / 2.0);
      double y = Sym(cavDim.y / 2.0);
      double z = Sym(cavDim.z / 2.0);
      loc.Set(i, x, y, z);
    }
    // Shift by center
    loc.AddRange(0, len, cavCenter);
  }

  // using UniformRandom algorithm in TransformMatrix.h
  XYZ RandomUnitVect() {
    double u2 = gen->rand(2.0 * M_PI);
    double u3 = gen->rand(2.0);
    double r = std::sqrt(u3);
    double root = std::sqrt(2.0 - u3);
    XYZ temp(std::sin(u2) * r * root, std::cos(u2) * r * root, 1.0 - u3);
    return temp;
  }

  void FillWithRandomOnSphere(XYZArray &loc, const uint len,
                              const double rAttach, const XYZ &center) {
    // Pick on cos(phi) - this was faster and always uses 2 rand calls
    for (uint i = 0; i < len; ++i) {
      XYZ sphere_val = PickOnUnitSphere();
      loc.Set(i, sphere_val * rAttach + center);
    }
  }

  // Returns a uniformly random point on the unit sphere
  XYZ PickOnUnitSphere() {
    // picking phi uniformly will cluster points at poles
    // pick u = cos(phi) uniformly instead
    double u = Sym(1.0);
    double theta = gen->randExc(2.0 * M_PI);
    double rootTerm = std::sqrt(1.0 - u * u);
    return XYZ(rootTerm * std::cos(theta), rootTerm * std::sin(theta), u);
  }

  // Pick using array of prob. -- used to pick move, weighted selection, etc.
  void PickArbDist(uint &pick, double &subDraw, double const *const w,
                   const double Wt, const uint len) {
    double sum = 0.0, prevSum = 0.0, draw = rand(Wt);
    for (pick = 0; pick < len && sum < draw; ++pick) {
      prevSum = sum;
      sum += w[pick];
    }
    if (pick != 0)
      --pick;
    subDraw = draw - prevSum;
  }

  // Pick an integer in range 0 -- n given a list of weights, and their sum
  // totalWeight
  uint PickWeighted(const double *weights, const int n,
                    const double totalWeight) {
    double draw = rand(totalWeight);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      sum += weights[i];
      if (sum >= draw) {
        return i;
      }
    }
    int lastZero = n;
    for (int i = n - 1; i > 0 && weights[i] == 0.0; --i) {
      lastZero = i;
    }
    lastZero--;
    if (std::fabs(draw - totalWeight) < 0.001) {
      return lastZero;
    }

    // If the code gets to here that means the total of all the weights was
    // more than totalWeight. So let's print out a message and exit the program
    std::cerr << "Error: In PRNG::PickWeighted() the total of all weights is\n";
    std::cerr << "more than totalWeight\nDebug info:\n";
    std::cerr << "totalWeight: " << totalWeight << "\n";
    std::cerr << "draw: " << draw << "\n";
    std::cerr << "sum: " << sum << "\n";

#if GOMC_LIB_MPI
    std::cerr
        << "Error: When conducting replica exchange simulations, this error\n";
    std::cerr
        << "usually occurs when using Intraswap or Molecular Exchange moves\n";
    std::cerr << "and seeding all simulations with the same seed.\n";
#endif

    exit(EXIT_FAILURE);
  }

  void PickBox(uint &b, const double subDraw, const double movPerc) const {
    // Calculate "chunk" of move for each box.
    double boxDiv = movPerc / BOX_TOTAL;
    // Which chunk was our draw in...
    b = (uint)(subDraw / boxDiv);
    // FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b >= BOX_TOTAL)
      b--;
  }

  void PickBox(uint &b, double &subDraw, double &boxDiv,
               const double movPerc) const {
    // Calculate "chunk" of move for each box.
    boxDiv = movPerc / BOX_TOTAL;
    // Which chunk was our draw in...
    b = (uint)(subDraw / boxDiv);
    // clamp if some rand err.
    // FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b >= BOX_TOTAL)
      b--;
    // Get draw down to box we picked, then calculate the size of
    // each molecule kind's chunk.
    subDraw -= b * boxDiv;
  }

  void PickBool(bool &b, const double subDraw, const double movPerc) const {
    // Calculate "chunk" of move for each box.
    double boxDiv = movPerc / 2;
    // Which chunk was our draw in...
    uint bpick = (uint)(subDraw / boxDiv);
    // FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (bpick >= BOX_TOTAL)
      bpick--;
    b = (bpick > 0 ? true : false);
  }

  void SetOtherBox(uint &bDest, const uint bSrc) const {
    // NOTE: assumes only two boxes.
    if (bSrc == mv::BOX0)
      bDest = mv::BOX1;
    else
      bDest = mv::BOX0;
  }

  // Function override (1 of 2)
  // Pick just a box pair.
  void PickBoxPair(uint &bSrc, uint &bDest, const double subDraw,
                   const double movPerc) {
    PickBox(bSrc, subDraw, movPerc);
    SetOtherBox(bDest, bSrc);
  }

  // Function override (2 of 2)
  // Pick a box pair, plus save leftovers to calculate mol. kind
  void PickBoxPair(uint &bSrc, uint &bDest, double &boxDiv, double &subDraw,
                   const double movPerc) {
    PickBox(bSrc, subDraw, boxDiv, movPerc);
    SetOtherBox(bDest, bSrc);
  }

  // Returns false if none of that kind of molecule in selected box or molecule
  // is fixed using tag (beta == 1).
  uint PickMol(uint &m, uint &mk, const uint b, const double subDraw,
               const double subPerc) {
    uint rejectState = mv::fail_state::NO_FAIL;
    uint mkTot = molLookRef.GetNumCanMoveKind();
    if (mkTot == 0) {
      std::cerr << "Error: All molecules inside the box are fixed!\n";
      exit(EXIT_FAILURE);
    }
    double molDiv = subPerc / mkTot;
    // Which molecule kind chunk are we in?
    uint k = (uint)(subDraw / molDiv);

    // clamp if some rand err.
    if (k == mkTot)
      k = mkTot - 1;

    mk = molLookRef.GetCanMoveKind(k);

    // Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      // Among the ones of that kind in that box, pick one @ random.
      // Molecule with a tag (beta == 1) cannot be selected.
      do {
        uint mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        // Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while (molLookRef.IsFix(m));
    }

    return rejectState;
  }

  // Returns false if none of that kind of molecule in selected box or molecule
  // is fixed using tag (beta == 1).
  uint PickMol(uint &m, uint &mk, const uint b) {
    uint rejectState = mv::fail_state::NO_FAIL;
    uint mkTot = molLookRef.GetNumCanMoveKind();
    if (mkTot == 0) {
      std::cerr << "Error: All molecules inside the box are fixed!\n";
      exit(EXIT_FAILURE);
    }
    // Which molecule kind chunk are we in?
    uint k = randIntExc(mkTot);

    // clamp if some rand err.
    if (k == mkTot)
      k = mkTot - 1;

    mk = molLookRef.GetCanMoveKind(k);

    // Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      // Among the ones of that kind in that box, pick one @ random.
      // Molecule with a tag (beta == 1) cannot be selected.
      do {
        uint mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        // Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while (molLookRef.IsFix(m));
    }

    return rejectState;
  }

  // Returns false if none of that kind of molecule in selected box or molecule
  // is fixed using tag (beta >= 1).
  uint PickMol2(uint &m, uint &mk, const uint b, const double subDraw,
                const double subPerc) {
    uint rejectState = mv::fail_state::NO_FAIL;
    uint mkTot = molLookRef.GetNumCanSwapKind();
    if (mkTot == 0) {
      std::cerr << "Error: All molecules inside the box are fixed!\n";
      exit(EXIT_FAILURE);
    }
    double molDiv = subPerc / mkTot;
    // Which molecule kind chunk are we in?
    uint k = (uint)(subDraw / molDiv);

    // clamp if some rand err.
    if (k == mkTot)
      k = mkTot - 1;

    mk = molLookRef.GetCanSwapKind(k);

    // Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      // Among the ones of that kind in that box, pick one @ random.
      // Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
      do {
        uint mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        // Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while (molLookRef.IsNoSwap(m));
    }

    return rejectState;
  }

  // Pick a random molecule of kind mk in box b
  uint PickMolIndex(uint &m_idx, const uint mk, const uint b) {
    uint rejectState = mv::fail_state::NO_FAIL;
    // Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      uint mOff, m;
      // Among the ones of that kind in that box, pick one @ random.
      // Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
      do {
        mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        // Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while (molLookRef.IsNoSwap(m));
      m_idx = m;
    }
    return rejectState;
  }

  // In MEMC move pick n molecules of kind mk
  uint PickMol(const uint mk, uint &mk2, uint &m2, const uint b) {
    uint rejectState = mv::fail_state::NO_FAIL;
    // Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      uint mOff, m;
      // Among the ones of that kind in that box, pick one @ random.
      // Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
      do {
        mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        // Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while (molLookRef.IsNoSwap(m));
      m2 = m;
      mk2 = mk;
    }

    return rejectState;
  }

  // In MEMC move pick n molecules of kind mk
  uint PickMol(const uint mk, std::vector<uint> &mk2, std::vector<uint> &m2,
               const uint n, const uint b) {
    uint rejectState = mv::fail_state::NO_FAIL;
    // Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else if (molLookRef.NumKindInBox(mk, b) < n) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      for (uint i = 0; i < n; i++) {
        uint mOff, m;
        // Among the ones of that kind in that box, pick one @ random.
        // Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
        do {
          mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
          // Lookup true index in table.
          m = molLookRef.GetMolNum(mOff, mk, b);
        } while (molLookRef.IsNoSwap(m) ||
                 std::find(m2.begin(), m2.end(), m) != m2.end());
        m2.push_back(m);
        mk2.push_back(mk);
      }
    }

    return rejectState;
  }

  // pick a molecule that is not fixed (beta != 1)
  uint PickMolAndBoxPair(uint &m, uint &mk, uint &bSrc, uint &bDest,
                         double subDraw, const double movPerc) {
    double boxDiv = 0;
    PickBoxPair(bSrc, bDest, boxDiv, subDraw, movPerc);
    return PickMol(m, mk, bSrc, subDraw, boxDiv);
  }

  // pick a molecule that can be transfer to other box (beta == 0)
  uint PickMolAndBoxPair2(uint &m, uint &mk, uint &bSrc, uint &bDest,
                          double subDraw, const double movPerc) {
    double boxDiv = 0;
    PickBoxPair(bSrc, bDest, boxDiv, subDraw, movPerc);
    return PickMol2(m, mk, bSrc, subDraw, boxDiv);
  }

  uint PickMolAndBox(uint &m, uint &mk, uint &b, double subDraw,
                     const double movPerc) {
    double boxDiv = 0;
    PickBox(b, subDraw, boxDiv, movPerc);
    return PickMol(m, mk, b, subDraw, boxDiv);
  }

  // draw a uniform distribution with average mean
  // and standard deviation stdDev
  double Gaussian(double mean, double stdDev) {
    if (hasSecondGaussian) {
      hasSecondGaussian = false;
      return (mean + secondGaussian * stdDev);
    } else {
      double r, v1, v2, factor;
      do {
        v1 = SymExc(1.0);
        v2 = SymExc(1.0);
        r = v1 * v1 + v2 * v2;
      } while (r >= 1.0 || r < 1.523e-8);

      factor = std::sqrt(-2.0 * std::log(r) / r);
      hasSecondGaussian = true;
      secondGaussian = v1 * factor;
      return (mean + v2 * factor * stdDev);
    }
  }

  MTRand *GetGenerator() { return gen; }

  bool operator==(const PRNG &other) { return (*gen == *other.gen); }

private:
  MTRand *gen;
  MoleculeLookup &molLookRef;
  bool hasSecondGaussian;
  double secondGaussian;
};

// Saves the current state of the PRNG as ./filename
inline void PRNG::saveState(const char *filename) {
  std::ofstream fout;
  fout.open(filename);

  if (!fout.is_open()) {
    std::cerr << "Failed to save PRNG state to file: " << filename << '\n';
    return;
  }

  MTRand::uint32 *saveArray = new MTRand::uint32[MTRand::N + 1];
  gen->save(saveArray);

  for (uint i = 0; i < MTRand::N + 1; ++i) {
    fout << saveArray[i] << '\n';
  }

  fout.close();

  delete[] saveArray;
}

#endif /*PRNG_H*/
