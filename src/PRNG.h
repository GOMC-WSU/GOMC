/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.0
Copyright (C) 2016  GOMC Group
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
    return (uint)(gen->randInt(bound-1));
  }

  //Generates number on (-bound,bound)
  double Sym(double bound)
  {
    return 2*bound*gen->rand() - bound;
  }

  /////////////////////////////
  //   GENERATION FUNCTIONS  //
  /////////////////////////////

  XYZ SymXYZ(double bound)
  {
    double bound2 = 2*bound;
    return XYZ(gen->rand(bound2)-bound, gen->rand(bound2)-bound,
               gen->rand(bound2)-bound);
  }

  //Used to pick first position of
  void FillWithRandom(XYZArray & loc, const uint len, XYZ const& axis)
  {
    for (uint i = 0; i < len; ++i)
      loc.Set(i, randExc(axis.x), randExc(axis.y), randExc(axis.z));
  }

  void FillWithRandomOnSphere(XYZArray & loc, const uint len,
                              const double rAttach, const XYZ& center)
  {
    //Quaternion Method - this was 80% slower in my tests - BGJ
    /*
    XYZ point;
    double x[4], sum;
    for (uint i = 0; i < len; ++i)
    {
    do
     {
        sum = 0;
        for (uint j = 0; j < 4; ++j)
        {
           x[j]=Sym(1.0);
           sum += x[j]*x[j];
        }
     } while (sum>=1);

     point.x = 2*(x[1]*x[3]+x[0]*x[2])/sum*rAttach;
     point.y = 2*(x[2]*x[3]-x[0]*x[1])/sum*rAttach;
     point.z = (x[0]*x[0]+x[3]*x[3]-x[1]*x[1]-x[2]*x[2])/sum*rAttach;
     loc.Set(i, point + center);
        }
        */
    //Pick on cos(phi) - this was faster and always uses 2 rand calls
    for (uint i = 0; i < len; ++i)
    {
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
    for (; pick < len && sum < draw; ++pick)
    {
      prevSum = sum;
      sum += w[pick];
    }
    if (pick!=0)
      --pick;
    subDraw = draw-prevSum;
  }

  //Pick an integer in range 0, n given a list of weights, and their sum totalWeight
  uint PickWeighted(const double *weights, const uint n, double totalWeight)
  {
    double draw = rand(totalWeight);
    double sum = 0.0;
    for(uint i = 0; i < n; ++i)
    {
      sum += weights[i];
      if(sum >= draw)
        return i;
    }
    //ya dun goofed, totalWeight was greater than the actual total
    //enjoy your segfault (this shouldn't fail gracefully)
    return -1;
  }


  void PickBox(uint &b,
               const double subDraw, const double movPerc) const
  {
    //Calculate "chunk" of move for each box.
    double boxDiv = movPerc/BOX_TOTAL;
    //Which chunk was our draw in...
    b = (uint)(subDraw/boxDiv);
    //FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b != mv::BOX0)
      b = mv::BOX1;
  }

  void PickBox(uint &b, double &subDraw, double &boxDiv,
               const double movPerc) const
  {
    //Calculate "chunk" of move for each box.
    boxDiv = movPerc/BOX_TOTAL;
    //Which chunk was our draw in...
    b = (uint)(subDraw/boxDiv);
    //clamp if some rand err.
    //FIXME: Hack to prevent picking invalid boxes, may violate balance
    if (b != mv::BOX0)
      b = mv::BOX1;
    //Get draw down to box we picked, then calculate the size of
    //each molecule kind's chunk.
    subDraw-= b*boxDiv;
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

  //Returns false if none of that kind of molecule in selected box.
  uint PickMol(uint & m, uint & mk, const uint b,
               const double subDraw, const double subPerc)
  {
    uint rejectState = mv::fail_state::NO_FAIL;
    uint mkTot = molLookRef.GetNumKind();
    double molDiv = subPerc/mkTot;
    //Which molecule kind chunk are we in?
    mk = (uint)(subDraw/molDiv);

    //clamp if some rand err.
    if (mk == mkTot)
      mk = mkTot -1;

    //Pick molecule with the help of molecule lookup table.
    if (molLookRef.NumKindInBox(mk, b) == 0)
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    else
    {
      //Among the ones of that kind in that box, pick one @ random.
      uint mOff = randIntExc(molLookRef.NumKindInBox(mk, b));
      //Lookup true index in table.
      m = molLookRef.GetMolNum(mOff, mk, b);
    }
    return rejectState;
  }

  uint PickMolAndBoxPair(uint &m, uint &mk, uint & bSrc, uint & bDest,
                         double subDraw, const double movPerc)
  {
    double boxDiv=0;
    PickBoxPair(bSrc, bDest, boxDiv, subDraw, movPerc);
    return PickMol(m, mk, bSrc, subDraw, boxDiv);
  }

  uint PickMolAndBox(uint & m, uint &mk, uint &b,
                     double subDraw, const double movPerc)
  {
    double boxDiv=0;
    PickBox(b, subDraw, boxDiv, movPerc);
    return PickMol(m, mk, b, subDraw, boxDiv);
  }

private:
  MTRand * gen;
  MoleculeLookup & molLookRef;
};

//Saves the current state of the PRNG as ./filename
inline void PRNG::saveState(const char* filename)
{
  std::ofstream fout;
  fout.open(filename);

  if(!fout.is_open())
  {
    std::cerr << "Failed to save PRNG state to file: " << filename << '\n';
    return;
  }

  MTRand::uint32* saveArray = new MTRand::uint32[MTRand::N + 1];
  gen->save(saveArray);

  for(uint i = 0; i < MTRand::N + 1; ++i)
  {
    fout << saveArray[i] << '\n';
  }

  fout.close();

  delete[] saveArray;

}

#endif /*PRNG_H*/
