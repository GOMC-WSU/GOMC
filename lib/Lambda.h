/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef LAMBDA_H
#define LAMBDA_H

//Defining lambda class to handle fractional molecule
class Lambda
{
public:
  Lambda() {
    std::fill_n(isFraction, BOX_TOTAL, false);
  }

  void Set(const double vdw, const double coulomb, const uint mol,
           const uint kind, const uint box);

  void UnSet(const uint sourceBox, const uint destBox);

  double GetLambdaVDW(const uint mol, const uint kind, const uint box) const;

  double GetLambdaVDW(const uint kind, const uint box) const;

  double GetLambdaCoulomb(const uint mol,const uint kind, const uint box) const;

  bool KindIsFractional(const uint kind, const uint box) const;

protected:
  uint molIndex[BOX_TOTAL];
  uint kindIndex[BOX_TOTAL];
  double lambdaVDW[BOX_TOTAL], lambdaCoulomb[BOX_TOTAL];
  bool isFraction[BOX_TOTAL];
};

inline void Lambda::Set(const double vdw, const double coulomb, const uint mol,
                const uint kind, const uint box)
{
    molIndex[box] = mol;
    kindIndex[box] = kind;
    lambdaVDW[box] = vdw;
    lambdaCoulomb[box] = coulomb;
    isFraction[box] = true;
}

inline void Lambda::UnSet(const uint sourceBox, const uint destBox) 
{
    isFraction[sourceBox] = isFraction[destBox] = false;
    lambdaVDW[sourceBox] = lambdaCoulomb[sourceBox] = 1.0;
    lambdaVDW[destBox] = lambdaCoulomb[destBox] = 0.0;
}

inline double Lambda::GetLambdaVDW(const uint mol, const uint kind,
                            const uint box) const
{
    double val = 1.0;
    if(isFraction[box]) {
        if((molIndex[box] == mol) && (kindIndex[box] == kind)) {
            val = lambdaVDW[box];
        }
    }
    return val;
}

inline double Lambda::GetLambdaVDW(const uint kind, const uint box) const
{
    double val = 1.0;
    if(isFraction[box]) {
        if(kindIndex[box] == kind) {
            val = lambdaVDW[box];
        }
    }
    return val;
}

inline double Lambda::GetLambdaCoulomb(const uint mol, const uint kind,
                                const uint box) const
{
    double val = 1.0;
    if(isFraction[box]) {
        if((molIndex[box] == mol) && (kindIndex[box] == kind)) {
            val = lambdaCoulomb[box];
        }
    }
    return val;
}

inline bool Lambda::KindIsFractional(const uint kind, const uint box) const 
{
    bool result = false;
    if(isFraction[box]) {
        if(kindIndex[box] == kind) {
            result = true;
        }
    }
    return result;
}

#endif