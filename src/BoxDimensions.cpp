/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BoxDimensions.h"
#include "MoveConst.h" //For cutoff-related fail condition

void BoxDimensions::Init(config_setup::RestartSettings const& restart,
                         config_setup::Volume const& confVolume,
                         pdb_setup::Cryst1 const& cryst,
                         double rc, double rcSq)
{
  const double TENTH_ANGSTROM = 0.1;
  rCut = rc;
  rCutSq = rcSq;
  minBoxSize = rc * rcSq * 8 + TENTH_ANGSTROM;
  if (restart.enable && cryst.hasVolume)
    axis = cryst.axis;
  else if (confVolume.hasVolume)
    axis = confVolume.axis;
  else
  {
    fprintf(stderr,
            "Error: Box Volume(s) not specified in PDB or in.dat files.\n");
    exit(EXIT_FAILURE);
  }
  halfAx.Init(BOX_TOTAL);
  axis.CopyRange(halfAx, 0, 0, BOX_TOTAL);
  halfAx.ScaleRange(0, BOX_TOTAL, 0.5);
  //Init volume/inverse volume.
  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    volume[b] = axis.x[b] * axis.y[b] * axis.z[b];
    volInv[b] = 1.0 / volume[b];
    //check to see if initial box size is cubic or not
    cubic[b] = ((axis.x[b] == axis.y[b]) && (axis.y[b] == axis.z[b]));
  }
  constArea = confVolume.cstArea;
}

uint BoxDimensions::ShiftVolume
(BoxDimensions & newDim, XYZ & scale, const uint b, const double delta) const
{
  uint rejectState = mv::fail_state::NO_FAIL;
  double newVolume = volume[b] + delta;

  newDim.SetVolume(b, newVolume);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  if ((newDim.axis.x[b] < rCut || newDim.axis.y[b] < rCut || newDim.axis.z[b] < rCut))
  {
    std::cout << "WARNING!!! box shrunk below 2*rc! Auto-rejecting!"
	      << std::endl;
    rejectState = mv::fail_state::VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF;
  }
  scale = newDim.axis.Get(b) / axis.Get(b);

  return rejectState;
}

uint BoxDimensions::ExchangeVolume
(BoxDimensions & newDim, XYZ * scale, const double transfer) const
{
  uint state = mv::fail_state::NO_FAIL;
  double vTot = volume[0] + volume[1];
  newDim = *this;

  newDim.SetVolume(0, volume[0] + transfer);
  newDim.SetVolume(1, vTot - newDim.volume[0]);

  //If move would shrink any box axis to be less than 2 * rcut, then
  //automatically reject to prevent errors.
  for (uint b = 0; b < BOX_TOTAL && state == mv::fail_state::NO_FAIL; b++)
  {
    scale[b] = newDim.axis.Get(b) / axis.Get(b);
    if ((newDim.halfAx.x[b] < rCut || newDim.halfAx.y[b] < rCut ||
	 newDim.halfAx.z[b] < rCut))
    {
      std::cout << "WARNING!!! box shrunk below 2*Rcut! Auto-rejecting!"
	      << std::endl;
      state = state && mv::fail_state::VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF;
    }
  }
  return state;
}
