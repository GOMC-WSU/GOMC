/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EWALDCACHED_H
#define EWALDCACHED_H

#include "Ewald.h"

class EwaldCached : public Ewald
{
public:

  EwaldCached(StaticVals & stat, System & sys);
  ~EwaldCached();

  virtual void Init();

  virtual void AllocMem();

  //setup reciprocate term for a box
  virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

  //calculate reciprocate energy term for a box
  virtual double BoxReciprocal(uint box) const;

  //calculate reciprocate term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
                               const uint box);

  //calculate reciprocate term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int molIndex);

  //calculate reciprocate term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                 const uint box, const int molIndex);
    
  //calculate reciprocate term for inserting some molecules (kindA) in
  //destination box and removing a molecule (kindB) from destination box
  virtual double SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
                           const std::vector<cbmc::TrialMol> &oldMol);

  //restore cosMol and sinMol
  virtual void RestoreMol(int molIndex);

  //update sinMol and cosMol
  virtual void exgMolCache();

  //backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache();

private:

  double *cosMolRestore; //cos()*charge
  double *sinMolRestore; //sin()*charge
  double **cosMolRef;
  double **sinMolRef;
  double **cosMolBoxRecip;
  double **sinMolBoxRecip;
  #if ENSEMBLE == GEMC
  const uint GEMC_KIND;
  #endif
};


#endif /*EWALDCACHED_H*/
