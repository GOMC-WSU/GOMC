
/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PRNG_H
#define PRNG_H

#include <fstream>
#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>

#include "MersenneTwister.h"
#include "EnsemblePreprocessor.h"
#include "MoleculeLookup.h"
#include "MoveConst.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "XYZArray.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "BasicTypes.h"


//Wrapper class for our random numbers
class PRNG
{
public:
  PRNG(MoleculeLookup & molLook) : molLookRef(molLook) {}
  ~PRNG(void)
  {
    delete gen;
  }

  void Init(MTRand * prng)
  {
    gen = prng;
  }

  //Saves the current state of the PRNG as ./filename
  void saveState(const char* filename);

  ////
  // BASIC GENERATION
  ////

  //Standard double generation on [0,1.0]
  double operator()()
  {
    return (*gen)();
  }

  //Generate a double on a [0,b]
  double rand(double const bound)
  {
    return gen->rand(bound);
  }

  //Generate a double on a [0,b)
  double randExc(double const bound)
  {
    return gen->randExc(bound);
  }

  //Generate an unsigned int on [0,bound]
  uint randInt(const uint bound)
  {
    return (uint)(gen->randInt(bound));
  }

  //Generate an unsigned int on [0,bound)
  uint randIntExc(const uint bound)
  {
    return (uint)(gen->randInt(bound - 1));
  }

  //Generates number on (-bound,bound)
  double Sym(double bound)
  {
    return 2 * bound * gen->rand() - bound;
  }

  /////////////////////////////
  //   GENERATION FUNCTIONS  //
  /////////////////////////////

  XYZ SymXYZ(double bound)
  {
    double bound2 = 2 * bound;
    return XYZ(gen->rand(bound2) - bound, gen->rand(bound2) - bound,
               gen->rand(bound2) - bound);
  }
    
  // return between [-bound, bound]
  double SymExc(double bound)
  {
    return 2 * gen->rand(bound) - bound;
  }

  //Used to pick first position of
  void FillWithRandom(XYZArray & loc, const uint len, BoxDimensions const& dims,
                      const uint b)
  {
    for (uint i = 0; i < len; ++i) {
      loc.Set(i, randExc(dims.axis[b][0]), randExc(dims.axis[b][1]),
              randExc(dims.axis[b][2]));
      loc.Set(i, dims.TransformSlant(loc.Get(i), b));
    }
  }

  void FillWithRandom(XYZ & loc, BoxDimensions const& dims,
                      const uint b)
  {
    XYZ temp(randExc(dims.axis[b][0]), randExc(dims.axis[b][1]),
	     randExc(dims.axis[b][2]));
    loc = dims.TransformSlant(temp, b);
  }

  void FillWithRandomInCavity(XYZ &loc, XYZ const& cavDim)
  {
    XYZ temp(SymExc(cavDim.x/2.0), SymExc(cavDim.y/2.0), SymExc(cavDim.z/2.0));
    loc = temp;
  }

  //using UniformRandom algorithm in TransformMatrix.h
  XYZ RandomUnitVect()
  {
    double u2 = gen->rand();
    double u3 = gen->rand();
    u2 *= 2.0 * M_PI;
    u3 *= 2.0;
    double r = sqrt(u3);
    double root = sqrt(2.0 - u3);
    XYZ temp(sin(u2) * r * root, cos(u2) * r * root, 1.0 - u3);
    return  temp;
  }

  void FillWithRandomOnSphere(XYZArray & loc, const uint len,
                              const double rAttach, const XYZ& center)
  {
    //Pick on cos(phi) - this was faster and always uses 2 rand calls
    for (uint i = 0; i < len; ++i) {
      loc.Set(i, PickOnUnitSphere() * rAttach + center);
    }
  }

  //Returns a uniformly random point on the unit sphere
  XYZ PickOnUnitSphere()
  {
    //picking phi uniformly will cluster points at poles
    //pick u = cos(phi) uniformly instead
    double u = gen->rand(2.0);
    u -= 1.0;
    double theta = gen->randExc(2 * M_PI);
    double rootTerm = sqrt(1 - u * u);
    return XYZ(rootTerm * cos(theta), rootTerm * sin(theta), u);
  }


  //Pick using array of prob. -- used to pick move, weighted selection, etc.
  void PickArbDist(uint & pick, double & subDraw,
                   double const*const w, const double Wt, const uint len)
  {
    double sum = 0.0, prevSum = 0.0, draw = rand(Wt);
    pick = 0;
    for (; pick < len && sum < draw; ++pick) {
      prevSum = sum;
      sum += w[pick];
    }
    if (pick != 0)
      --pick;
    subDraw = draw - prevSum;
  }

  //Pick an integer in range 0, n given a list of weights, and their sum totalWeight
  uint PickWeighted(const double *weights, const uint n, double totalWeight)
  {
    double draw = rand(totalWeight);
    double sum = 0.0;
    for(uint i = 0; i < n; ++i) {
      sum += weights[i];
      if(sum >= draw) {
        return i;
      }
    }
    int lastZero = n;
    for(uint i = n - 1; i > 0 && !weights[i]; i--) {
      lastZero = i;
    }
    lastZero--;
    if(abs(draw - totalWeight) < 0.001) {
      return lastZero;
    }

    // If the code gets to here that means the total of all the weights was
    // more than totalWeight. So let's print out a message and exit the program
    std::cerr << "Error: In PRNG::PickWeighted() the total of all weights was" << std::endl;
    std::cerr << "more than totalWeight" << std::endl << "Debug info:\n";
    std::cerr << "totalWeight: " << totalWeight << std::endl;
    std::cerr << "draw: " << draw << std::endl;
    std::cerr << "sum: " << sum << std::endl;
    exit(EXIT_FAILURE);
  }


  void PickBox(uint &b, const double subDraw, const double movPerc) const
  {
    //Calculate "chunk" of move for each box.
    double boxDiv = movPerc / BOX_TOTAL;
    //Which chunk was our draw in...
    b = (uint)(subDraw / boxDiv);
    //FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b != mv::BOX0)
      b = mv::BOX1;
  }

  void PickBool(bool &b, const double subDraw, const double movPerc) const
  {
    //Calculate "chunk" of move for each box.
    double boxDiv = movPerc / 2;
    //Which chunk was our draw in...
    b = (uint)(subDraw/boxDiv);
    //FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b != mv::BOX0)
      b = mv::BOX1;
    b = (b > 0 ? true : false);
  }

  void PickBox(uint &b, double &subDraw, double &boxDiv,
               const double movPerc) const
  {
    //Calculate "chunk" of move for each box.
    boxDiv = movPerc / BOX_TOTAL;
    //Which chunk was our draw in...
    b = (uint)(subDraw / boxDiv);
    //clamp if some rand err.
    //FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b != mv::BOX0)
      b = mv::BOX1;
    //Get draw down to box we picked, then calculate the size of
    //each molecule kind's chunk.
    subDraw -= b * boxDiv;
  }

  void SetOtherBox(uint & bDest, const uint bSrc) const
  {
    //NOTE: assumes only two boxes.
    if (bSrc == mv::BOX0)
      bDest = mv::BOX1;
    else
      bDest = mv::BOX0;
  }

  //Function override (1 of 2)
  //Pick just a box pair.
  void PickBoxPair(uint & bSrc, uint & bDest,
                   const double subDraw, const double movPerc)
  {
    PickBox(bSrc, subDraw, movPerc);
    SetOtherBox(bDest, bSrc);
  }

  //Function override (2 of 2)
  //Pick a box pair, plus save leftovers to calculate mol. kind
  void PickBoxPair(uint & bSrc, uint & bDest, double & boxDiv,
                   double & subDraw, const double movPerc)
  {
    PickBox(bSrc, subDraw, boxDiv, movPerc);
    SetOtherBox(bDest, bSrc);
  }

  //Returns false if none of that kind of molecule in selected box or molecule
  //is fixed using tag (beta == 1).
  uint PickMol(uint & m, uint & mk, const uint b,
               const double subDraw, const double subPerc)
  {
    uint rejectState = mv::fail_state::NO_FAIL;
    uint mkTot = molLookRef.GetNumCanMoveKind();
    double molDiv = subPerc / mkTot;
    //Which molecule kind chunk are we in?
    uint k = (uint)(subDraw / molDiv);

    //clamp if some rand err.
    if (k == mkTot)
      k = mkTot - 1;

    mk = molLookRef.GetCanMoveKind(k);

    //Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      //Among the ones of that kind in that box, pick one @ random.
      //Molecule with a tag (beta == 1) cannot be selected.
      do {
        uint mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        //Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while(molLookRef.IsFix(m));
    }

    return rejectState;
  }

  //Returns false if none of that kind of molecule in selected box or molecule
  //is fixed using tag (beta >= 1).
  uint PickMol2(uint & m, uint & mk, const uint b,
                const double subDraw, const double subPerc)
  {
    uint rejectState = mv::fail_state::NO_FAIL;
    uint mkTot = molLookRef.GetNumCanSwapKind();
    double molDiv = subPerc / mkTot;
    //Which molecule kind chunk are we in?
    uint k = (uint)(subDraw / molDiv);

    //clamp if some rand err.
    if (k == mkTot)
      k = mkTot - 1;

    mk = molLookRef.GetCanSwapKind(k);

    //Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      //Among the ones of that kind in that box, pick one @ random.
      //Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
      do {
        uint mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
        //Lookup true index in table.
        m = molLookRef.GetMolNum(mOff, mk, b);
      } while(molLookRef.IsNoSwap(m));
    }

    return rejectState;
  }

  // In MEMC move pick n molecules of kind mk
  uint PickMol(const uint mk, uint &mk2, uint &m2, const uint b)
  {
    uint rejectState = mv::fail_state::NO_FAIL;
    //Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      uint mOff, m;
      //Among the ones of that kind in that box, pick one @ random.
      //Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
      do {
	mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
	//Lookup true index in table.
	m = molLookRef.GetMolNum(mOff, mk, b);
      } while(molLookRef.IsNoSwap(m));
      m2 = m;
      mk2 = mk;
    }
    
    return rejectState;
  }

  // In MEMC move pick n molecules of kind mk
  uint PickMol(const uint mk, std::vector<uint> &mk2,
               std::vector<uint> &m2, const uint n, const uint b)
  {
    uint rejectState = mv::fail_state::NO_FAIL;
    //Pick molecule with the help of molecule lookup table.
    if ((molLookRef.NumKindInBox(mk, b) == 0)) {
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else if (molLookRef.NumKindInBox(mk, b) < n){
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    } else {
      for(uint i = 0; i < n; i++) {
        uint mOff, m;
        //Among the ones of that kind in that box, pick one @ random.
        //Molecule with a tag (beta == 2 and beta == 1) cannot be selected.
        do {
          mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
          //Lookup true index in table.
          m = molLookRef.GetMolNum(mOff, mk, b);
        } while(molLookRef.IsNoSwap(m) || 
		std::find(m2.begin(), m2.end(), m) != m2.end());
        m2.push_back(m);
        mk2.push_back(mk);
      }
    }
      
    return rejectState;
  }
    
  // pick a molecule that is not fixed (beta != 1)
  uint PickMolAndBoxPair(uint &m, uint &mk, uint & bSrc, uint & bDest,
                         double subDraw, const double movPerc)
  {
    double boxDiv = 0;
    PickBoxPair(bSrc, bDest, boxDiv, subDraw, movPerc);
    return PickMol(m, mk, bSrc, subDraw, boxDiv);
  }

  // pick a molecule that can be transfer to other box (beta == 0)
  uint PickMolAndBoxPair2(uint &m, uint &mk, uint & bSrc, uint & bDest,
                          double subDraw, const double movPerc)
  {
    double boxDiv = 0;
    PickBoxPair(bSrc, bDest, boxDiv, subDraw, movPerc);
    return PickMol2(m, mk, bSrc, subDraw, boxDiv);
  }

  uint PickMolAndBox(uint & m, uint &mk, uint &b,
                     double subDraw, const double movPerc)
  {
    double boxDiv = 0;
    PickBox(b, subDraw, boxDiv, movPerc);
    return PickMol(m, mk, b, subDraw, boxDiv);
  }

  MTRand * GetGenerator() { return gen; }

private:
  MTRand * gen;
  MoleculeLookup & molLookRef;
};

//Saves the current state of the PRNG as ./filename
inline void PRNG::saveState(const char* filename)
{
  std::ofstream fout;
  fout.open(filename);

  if(!fout.is_open()) {
    std::cerr << "Failed to save PRNG state to file: " << filename << '\n';
    return;
  }

  MTRand::uint32* saveArray = new MTRand::uint32[MTRand::N + 1];
  gen->save(saveArray);

  for(uint i = 0; i < MTRand::N + 1; ++i) {
    fout << saveArray[i] << '\n';
  }

  fout.close();

  delete[] saveArray;

}

#endif /*PRNG_H*/
