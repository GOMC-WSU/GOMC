/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <map> //for function handle storage.
#include <string> //for var names, etc.
#include <vector>
#include <string>
#include <iomanip>

#include "ConfigSetup.h"

//#define UINT_MAX 0xffffffff
//#define ULONG_MAX 0xffffffffUL
#define DBL_MAX 1.7976931348623158e+308

int stringtoi(const std::string& s)
{
  std::istringstream str(s);
  int i;
  str >> i;
  return i;
}

double stringtod(const std::string& s)
{
  std::istringstream str(s);
  double i;
  str >> i;
  return i;
}


ConfigSetup::ConfigSetup(void)
{
  int i;
  in.restart.enable = false;
  in.restart.step = ULONG_MAX;
  in.prng.seed = UINT_MAX;
  sys.elect.readEwald = false;
  sys.elect.readElect = false;
  sys.elect.readCache = false;
  sys.elect.ewald = false;
  sys.elect.enable = false;
  sys.elect.tolerance = DBL_MAX;
  sys.elect.oneFourScale = DBL_MAX;
  sys.elect.dielectric = DBL_MAX;
  sys.step.total = ULONG_MAX;
  sys.step.equil = ULONG_MAX;
  sys.step.adjustment = ULONG_MAX;
  sys.step.pressureCalcFreq = ULONG_MAX;
  sys.step.pressureCalc = true;
  in.ffKind.numOfKinds = 0;
  sys.exclude.EXCLUDE_KIND = UINT_MAX;
  in.prng.kind = "";
  in.files.param.name = "";
  for(i = 0; i < BOX_TOTAL; i++) {
    in.files.pdb.name[i] = "";
    in.files.psf.name[i] = "";
  }
#if ENSEMBLE == GEMC
  sys.gemc.kind = UINT_MAX;
  sys.gemc.pressure = DBL_MAX;
#endif
#if ENSEMBLE == NPT
  sys.gemc.kind = mv::GEMC_NPT;
  sys.gemc.pressure = DBL_MAX;
#endif

  sys.T.inKelvin = DBL_MAX;
  sys.ff.VDW_KIND = UINT_MAX;
  sys.ff.doTailCorr = true;
  sys.ff.rswitch = DBL_MAX;
  sys.ff.cutoff = DBL_MAX;
  sys.ff.cutoffLow = DBL_MAX;
  sys.moves.displace = DBL_MAX;
  sys.moves.rotate = DBL_MAX;
  sys.moves.intraSwap = DBL_MAX;
  sys.moves.multiParticle = DBL_MAX;
  out.state.settings.enable = true;
  out.restart.settings.enable = true;
  out.console.enable = true;
  out.statistics.settings.block.enable = true;
#if ENSEMBLE == GCMC
  sys.chemPot.isFugacity = false;
  out.statistics.settings.hist.enable = false;
  out.statistics.settings.hist.frequency = ULONG_MAX;
  out.state.files.hist.histName = "";
  out.state.files.hist.letter = "";
  out.state.files.hist.number = "";
  out.state.files.hist.sampleName = "";
  out.state.files.hist.stepsPerHistSample = UINT_MAX;
#endif
  out.statistics.settings.uniqueStr.val = "";
  out.state.settings.frequency = ULONG_MAX;
  out.restart.settings.frequency = ULONG_MAX;
  out.console.frequency = ULONG_MAX;
  out.statistics.settings.block.frequency = ULONG_MAX;
  out.statistics.vars.energy.block = false;
  out.statistics.vars.energy.fluct = false;
  out.statistics.vars.pressure.block = false;
  out.statistics.vars.pressure.fluct = false;
  out.statistics.vars.surfaceTension.block = false;
  out.statistics.vars.surfaceTension.fluct = false;
#ifdef VARIABLE_PARTICLE_NUMBER
  sys.moves.transfer = DBL_MAX;
  sys.cbmcTrials.bonded.ang = UINT_MAX;
  sys.cbmcTrials.bonded.dih = UINT_MAX;
  sys.cbmcTrials.nonbonded.first = UINT_MAX;
  sys.cbmcTrials.nonbonded.nth = UINT_MAX;
  out.statistics.vars.molNum.block = false;
  out.statistics.vars.molNum.fluct = false;
#endif
#ifdef VARIABLE_VOLUME
  sys.moves.volume = DBL_MAX;
  out.statistics.vars.volume.block = false;
  out.statistics.vars.volume.fluct = false;
#endif
  out.statistics.vars.density.block = false;
  out.statistics.vars.density.fluct = false;
}

bool ConfigSetup::checkBool(string str)
{
  int k;
  // capitalize string
  for(k = 0; k < str.length(); k++) {
    str[k] = toupper(str[k]);
  }

  if(str == "ON" || str == "TRUE" || str == "YES")
    return true;
  else if(str == "OFF" || str == "FALSE" || str == "NO")
    return false;
  std::cout << "Error: " << str << "couldn't be recognized!" << std::endl;
  exit(EXIT_FAILURE);
}

void ConfigSetup::Init(const char *fileName)
{
  std::vector<std::string> line;

  reader.Open(fileName);
  printf("\n%-40s %-s\n", "Reading Input File:", fileName);
  while(reader.readNextLine(line)) {
    if(line.size() == 0)
      continue;

    if(line[0] == "Restart") {
      in.restart.enable = checkBool(line[1]);
      if(in.restart.enable) {
        printf("%-40s %-s \n", "Info: Restart simulation",  "Active");
      }
    } else if(line[0] == "FirstStep") {
      in.restart.step = stringtoi(line[1]);
    } else if(line[0] == "PRNG") {
      in.prng.kind = line[1];
      if("RANDOM" == line[1])
        printf("%-40s %-s \n", "Info: Random seed", "Active");
    } else if(line[0] == "ParaTypeCHARMM") {
      if(checkBool(line[1])) {
        in.ffKind.numOfKinds++;
        in.ffKind.isEXOTIC = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isCHARMM = true;
        printf("Info: PARAMETER file: CHARMM format!\n");
      }
    } else if(line[0] == "ParaTypeEXOTIC") {
      if(checkBool(line[1])) {
        in.ffKind.numOfKinds++;
        in.ffKind.isCHARMM = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isEXOTIC = true;
        printf("Info: PARAMETER file: EXOTIC format!\n");
      }
    } else if(line[0] == "ParaTypeMARTINI") {
      if(checkBool(line[1])) {
        in.ffKind.numOfKinds ++;
        in.ffKind.isEXOTIC = false;
        in.ffKind.isMARTINI = true;
        in.ffKind.isCHARMM = true;
        printf("Info: PARAMETER file: MARTINI using CHARMM format!\n");
      }
    } else if(line[0] == "Parameters") {
      in.files.param.name = line[1];
    } else if(line[0] == "Coordinates") {
      uint boxnum = stringtoi(line[1]);
      if(boxnum >= BOX_TOTAL) {
        std::cout << "Error: Simulation requires " << BOX_TOTAL << " PDB file(s)!\n";
        exit(EXIT_FAILURE);
      }
      in.files.pdb.name[boxnum] = line[2];
    } else if(line[0] == "Structure") {
      uint boxnum = stringtoi(line[1]);
      if(boxnum >= BOX_TOTAL) {
        std::cout << "Error: Simulation requires " << BOX_TOTAL << " PSF file(s)!\n";
        exit(EXIT_FAILURE);
      }
      in.files.psf.name[boxnum] = line[2];
    }
#if ENSEMBLE == GEMC
    else if(line[0] == "GEMC") {
      if(line[1] == "NVT") {
        sys.gemc.kind = mv::GEMC_NVT;
        printf("Info: Running NVT_GEMC\n");
      } else if(line[1] == "NPT") {
        sys.gemc.kind = mv::GEMC_NPT;
        printf("Info: Running NPT_GEMC\n");
      }
    } else if(line[0] == "Pressure") {
      sys.gemc.pressure = stringtod(line[1]);
      printf("%-40s %-4.4f bar\n", "Info: Input Pressure", sys.gemc.pressure);
      sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
    }
#endif
#if ENSEMBLE == NPT
    else if(line[0] == "Pressure") {
      sys.gemc.kind = mv::GEMC_NPT;
      sys.gemc.pressure = stringtod(line[1]);
      printf("%-40s %-4.4f bar\n", "Info: Input Pressure", sys.gemc.pressure);
      sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
    }
#endif
    else if(line[0] == "Temperature") {
      sys.T.inKelvin = stringtod(line[1]);
      printf("%-40s %-4.4f K\n", "Info: Input Temperature", sys.T.inKelvin);
    } else if(line[0] == "Potential") {
      if(line[1] == "VDW") {
        sys.ff.VDW_KIND = sys.ff.VDW_STD_KIND;
        printf("%-40s %-s \n", "Info: Non-truncated potential", "Active");
      } else if(line[1] == "SHIFT") {
        sys.ff.VDW_KIND = sys.ff.VDW_SHIFT_KIND;
        printf("%-40s %-s \n", "Info: Shift truncated potential", "Active");
      } else if(line[1] == "SWITCH") {
        sys.ff.VDW_KIND = sys.ff.VDW_SWITCH_KIND;
        printf("%-40s %-s \n", "Info: Switch truncated potential", "Active");
      }
    } else if(line[0] == "LRC") {
      sys.ff.doTailCorr = checkBool(line[1]);
      printf("%-40s %-s \n", "Info: Long Range Correction", "Active");
    } else if(line[0] == "Rswitch") {
      sys.ff.rswitch = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Switch distance",
             sys.ff.rswitch);
    } else if(line[0] == "Rcut") {
      sys.ff.cutoff = stringtod(line[1]);
      printf("%-40s %-4.4f A\n", "Info: Cutoff", sys.ff.cutoff);
    } else if(line[0] == "RcutLow") {
      sys.ff.cutoffLow = stringtod(line[1]);
      printf("%-40s %-4.4f A.\n", "Info: Short Range Cutoff", sys.ff.cutoffLow);
    } else if(line[0] == "Exclude") {
      if(line[1] == sys.exclude.EXC_ONETWO) {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETWO_KIND;
        printf("%-40s %-s \n", "Info: Exclude", "ONE-TWO");
      } else if(line[1] == sys.exclude.EXC_ONETHREE) {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETHREE_KIND;
        printf("%-40s %-s \n", "Info: Exclude", "ONE-THREE");
      } else if(line[1] == sys.exclude.EXC_ONEFOUR) {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONEFOUR_KIND;
        printf("%-40s %-s \n", "Info: Exclude", "ONE-FOUR");
        printf("Warning: Modified 1-4 VDW parameters will be ignored!\n");
      }
    } else if(line[0] == "Ewald") {
      sys.elect.ewald = checkBool(line[1]);
      sys.elect.readEwald = true;
      if(sys.elect.ewald) {
        printf("%-40s %-s \n", "Info: Ewald Summation", "Active");
      }
    } else if(line[0] == "ElectroStatic") {
      sys.elect.enable = checkBool(line[1]);
      sys.elect.readElect = true;
    } else if(line[0] == "Tolerance") {
      sys.elect.tolerance = stringtod(line[1]);
      sys.elect.alpha = sqrt(-1 * log(sys.elect.tolerance)) /
                        sys.ff.cutoff;
      sys.elect.recip_rcut = 2 * (-log(sys.elect.tolerance)) /
                             sys.ff.cutoff;
      printf("%-40s %-1.3E \n", "Info: Ewald Summation Tolerance",
             sys.elect.tolerance);
    } else if(line[0] == "CachedFourier") {
      sys.elect.cache = checkBool(line[1]);
      sys.elect.readCache = true;
      if(sys.elect.cache) {
        printf("%-40s %-s \n", "Info: Cache Ewald Fourier", "Active");
      } else {
        printf("%-40s %-s \n", "Info: Cache Ewald Fourier", "Inactive");
      }
    } else if(line[0] == "1-4scaling") {
      sys.elect.oneFourScale = stringtod(line[1]);
    } else if(line[0] == "Dielectric") {
      sys.elect.dielectric = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Dielectric", sys.elect.dielectric);
    } else if(line[0] == "RunSteps") {
      sys.step.total = stringtoi(line[1]);
      printf("%-40s %-lu \n", "Info: Total number of steps", sys.step.total);
    } else if(line[0] == "EqSteps") {
      sys.step.equil = stringtoi(line[1]);
      printf("%-40s %-lu \n", "Info: Number of equilibration steps",
             sys.step.equil);
    } else if(line[0] == "AdjSteps") {
      sys.step.adjustment = stringtoi(line[1]);
      printf("%-40s %-lu \n", "Info: Move adjustment frequency",
             sys.step.adjustment);
    } else if(line[0] == "PressureCalc") {
      sys.step.pressureCalc = checkBool(line[1]);
      if(line.size() == 3)
        sys.step.pressureCalcFreq = stringtoi(line[2]);

      if(sys.step.pressureCalc && (line.size() == 2)) {
        std::cout <<
          "Error: Pressure calculation frequency is not specified!\n";
        exit(EXIT_FAILURE);
      }
      if(!sys.step.pressureCalc)
        printf("%-40s %-s \n", "Info: Pressure calculation", "Inactive");
      else {
        printf("%-40s %-lu \n", "Info: Pressure calculation frequency",
               sys.step.pressureCalcFreq);
      }
    } else if(line[0] == "DisFreq") {
      sys.moves.displace = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Displacement move frequency",
             sys.moves.displace);
    } else if(line[0] == "MultiParticleFreq") {
      sys.moves.multiParticle = stringtod(line[1]);
      printf("%-40s %-4.4f \n",
             "Info: Multi-Particle move frequency",
             sys.moves.multiParticle);
    } else if(line[0] == "IntraSwapFreq") {
      sys.moves.intraSwap = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Intra-Swap move frequency",
             sys.moves.intraSwap);
    } else if(line[0] == "RotFreq") {
      sys.moves.rotate = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Rotation move frequency",
             sys.moves.rotate);
    }
#ifdef VARIABLE_VOLUME
    else if(line[0] == "VolFreq") {
      sys.moves.volume = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Volume move frequency",
             sys.moves.volume);
    } else if(line[0] == "useConstantArea") {
      sys.volume.cstArea = checkBool(line[1]);
      if(sys.volume.cstArea)
        printf("Info: Volume change using constant X-Y area.\n");
      else
        printf("Info: Volume change using constant ratio.\n");
    } else if(line[0] == "FixVolBox0") {
      sys.volume.cstVolBox0 = checkBool(line[1]);
      if (sys.volume.cstVolBox0)
        printf("%-40s %-d \n", "Info: Fix volume box", 1);
    }
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
    else if(line[0] == "SwapFreq") {
#if ENSEMBLE == NVT || ENSEMBLE == NPT
      sys.moves.transfer = 0.000;
#else
      sys.moves.transfer = stringtod(line[1]);
#endif
      printf("%-40s %-4.4f \n", "Info: Molecule swap move frequency",
             sys.moves.transfer);
    }
#endif
    else if(line[0] == "CellBasisVector1") {
      uint box = stringtoi(line[1]);
      if(box < BOX_TOTAL) {
        XYZ temp;
        sys.volume.boxCoordRead++;
        sys.volume.hasVolume = (sys.volume.boxCoordRead == 3 * BOX_TOTAL);
        temp.x = stringtod(line[2]);
        temp.y = stringtod(line[3]);
        temp.z = stringtod(line[4]);
        sys.volume.axis[box].Set(0, temp);
        printf("%s %-d: %-26s %6.3f %7.3f %7.3f \n",
               "Info: Box ", box, " Periodic Cell Basis 1",
               temp.x, temp.y, temp.z);
      } else {
        std::cout << "Error: This simulation requires only " << BOX_TOTAL <<
                  " number of box dimension(s)!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else if(line[0] == "CellBasisVector2") {
      uint box = stringtoi(line[1]);
      if(box < BOX_TOTAL) {
        XYZ temp;
        sys.volume.boxCoordRead++;
        sys.volume.hasVolume = (sys.volume.boxCoordRead == 3 * BOX_TOTAL);
        temp.x = stringtod(line[2]);
        temp.y = stringtod(line[3]);
        temp.z = stringtod(line[4]);
        sys.volume.axis[box].Set(1, temp);
        printf("%s %-d: %-26s %6.3f %7.3f %7.3f \n",
               "Info: Box ", box, " Periodic Cell Basis 2",
               temp.x, temp.y, temp.z);
      } else {
        std::cout << "Error: This simulation requires only " << BOX_TOTAL <<
                  " number of box dimension(s)!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else if(line[0] == "CellBasisVector3") {
      uint box = stringtoi(line[1]);
      if(box < BOX_TOTAL) {
        XYZ temp;
        sys.volume.boxCoordRead++;
        sys.volume.hasVolume = (sys.volume.boxCoordRead == 3 * BOX_TOTAL);
        temp.x = stringtod(line[2]);
        temp.y = stringtod(line[3]);
        temp.z = stringtod(line[4]);
        sys.volume.axis[box].Set(2, temp);
        printf("%s %-d: %-26s %6.3f %7.3f %7.3f \n",
               "Info: Box ", box, " Periodic Cell Basis 3",
               temp.x, temp.y, temp.z);
      } else {
        std::cout << "Error: This simulation requires only " << BOX_TOTAL <<
                  " number of box dimension(s)!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#ifdef VARIABLE_PARTICLE_NUMBER
    else if(line[0] == "CBMC_First") {
      sys.cbmcTrials.nonbonded.first = stringtoi(line[1]);
    } else if(line[0] == "CBMC_Nth") {
      sys.cbmcTrials.nonbonded.nth = stringtoi(line[1]);
    } else if(line[0] == "CBMC_Ang") {
      sys.cbmcTrials.bonded.ang = stringtoi(line[1]);
    } else if(line[0] == "CBMC_Dih") {
      sys.cbmcTrials.bonded.dih = stringtoi(line[1]);
    }
#endif
#if ENSEMBLE == GCMC
    else if(line[0] == "ChemPot") {
      if(line.size() != 3) {
        std::cout << "Error: Chemical potential parameters are not specified!\n";
        exit(EXIT_FAILURE);
      }
      std::string resName = line[1];
      double val = stringtod(line[2]);
      sys.chemPot.cp[resName] = val;
      printf("%-40s %-6s %-6.4f K\n", "Info: Chemical potential",
             resName.c_str(), val);
    } else if(line[0] == "Fugacity") {
      if(line.size() != 3) {
        std::cout << "Error: Fugacity parameters are not specified!\n";
        exit(EXIT_FAILURE);
      }
      sys.chemPot.isFugacity = true;
      std::string resName = line[1];
      double val = stringtod(line[2]);
      sys.chemPot.cp[resName] = val * unit::BAR_TO_K_MOLECULE_PER_A3;
      printf("%-40s %-6s %-6.4f bar\n", "Info: Fugacity", resName.c_str(),
             val);
    }
#endif
    else if(line[0] == "OutputName") {
      out.statistics.settings.uniqueStr.val = line[1];
      printf("%-40s %-s \n", "Info: Output name", line[1].c_str());
    } else if(line[0] == "CoordinatesFreq") {
      out.state.settings.enable = checkBool(line[1]);
      if(line.size() == 3)
        out.state.settings.frequency = stringtoi(line[2]);

      if(out.state.settings.enable && (line.size() == 2))
        out.state.settings.frequency = (ulong)sys.step.total / 10;

      if(out.state.settings.enable) {
        printf("%-40s %-lu \n", "Info: Coordinate frequency",
               out.state.settings.frequency);
      } else
        printf("%-40s %-s \n", "Info: Printing coordinate", "Inactive");
    } else if(line[0] == "RestartFreq") {
      out.restart.settings.enable = checkBool(line[1]);
      if(line.size() == 3)
        out.restart.settings.frequency = stringtoi(line[2]);

      if(out.restart.settings.enable && (line.size() == 2))
        out.restart.settings.frequency = (ulong)sys.step.total;

      if(out.restart.settings.enable) {
        printf("%-40s %-lu \n", "Info: Restart frequency",
               out.restart.settings.frequency);
      } else
        printf("%-40s %-s \n", "Info: Printing restart coordinate", "Inactive");
    } else if(line[0] == "ConsoleFreq") {
      out.console.enable = checkBool(line[1]);
      if(line.size() == 3)
        out.console.frequency = stringtoi(line[2]);

      if(out.console.enable && (line.size() == 2)) {
        if(sys.step.total > 1000) {
          out.console.frequency = (ulong)sys.step.total / 1000;
        } else {
          out.console.frequency = (ulong)sys.step.total / 100;
        }
      }
      if(out.console.enable) {
        printf("%-40s %-lu \n", "Info: Console output frequency",
               out.console.frequency);
      } else
        printf("%-40s %-s \n", "Info: Console output", "Inactive");
    } else if(line[0] == "BlockAverageFreq") {
      out.statistics.settings.block.enable = checkBool(line[1]);
      if(line.size() == 3)
        out.statistics.settings.block.frequency = stringtoi(line[2]);

      if(out.statistics.settings.block.enable && (line.size() == 2))
        out.statistics.settings.block.frequency = (ulong)sys.step.total / 100;

      if(out.statistics.settings.block.enable) {
        printf("%-40s %-lu \n", "Info: Average output frequency",
               out.statistics.settings.block.frequency);
      } else
        printf("%-40s %-s \n", "Info: Average output", "Inactive");
    }
#if ENSEMBLE == GCMC
    else if(line[0] == "HistogramFreq") {
      out.statistics.settings.hist.enable = checkBool(line[1]);
      if(line.size() == 3)
        out.statistics.settings.hist.frequency = stringtoi(line[2]);

      if(out.statistics.settings.hist.enable && (line.size() == 2)) {
        if(sys.step.total > 1000) {
          out.statistics.settings.hist.frequency = (ulong)sys.step.total / 1000;
        } else {
          out.statistics.settings.hist.frequency = (ulong)sys.step.total / 100;
        }
      }

      if(out.statistics.settings.hist.enable) {
        printf("%-40s %-lu \n", "Info: Histogram output frequency",
               out.statistics.settings.hist.frequency);
      } else
        printf("%-40s %-s \n", "Info: Histogram output", "Inactive");
    } else if(line[0] == "DistName") {
      out.state.files.hist.histName = line[1];
    } else if(line[0] == "HistName") {
      out.state.files.hist.sampleName = line[1];
    } else if(line[0] == "RunNumber") {
      out.state.files.hist.number = line[1];
    } else if(line[0] == "RunLetter") {
      out.state.files.hist.letter = line[1];
    } else if(line[0] == "SampleFreq") {
      out.state.files.hist.stepsPerHistSample = stringtoi(line[1]);
      printf("%-40s %-d \n", "Info: Histogram sample frequency",
             out.state.files.hist.stepsPerHistSample);
    }
#endif
    else if(line[0] == "OutEnergy") {
      out.statistics.vars.energy.block = checkBool(line[1]);
      out.statistics.vars.energy.fluct = checkBool(line[2]);
    } else if(line[0] == "OutPressure") {
      out.statistics.vars.pressure.block = checkBool(line[1]);
      out.statistics.vars.pressure.fluct = checkBool(line[2]);
    }
#ifdef VARIABLE_PARTICLE_NUMBER
    else if(line[0] == "OutMolNum") {
      out.statistics.vars.molNum.block = checkBool(line[1]);
      out.statistics.vars.molNum.fluct = checkBool(line[2]);
    }
#endif
    else if(line[0] == "OutDensity") {
      out.statistics.vars.density.block = checkBool(line[1]);
      out.statistics.vars.density.fluct = checkBool(line[2]);
    } else if(line[0] == "OutSurfaceTension") {
      out.statistics.vars.surfaceTension.block = checkBool(line[1]);
      out.statistics.vars.surfaceTension.fluct = checkBool(line[2]);
    }
#ifdef VARIABLE_VOLUME
    else if(line[0] == "OutVolume") {
      out.statistics.vars.volume.block = checkBool(line[1]);
      out.statistics.vars.volume.fluct = checkBool(line[2]);
    }
#endif
    else if(line[0] == "Random_Seed") {
      in.prng.seed = stringtoi(line[1]);
      if("INTSEED" == in.prng.kind)
        printf("%-40s %-s \n", "Info: Constant seed", "Active");
      else
        printf("Warning: Constant seed set, but will be ignored.\n");
    } else {
      cout << "Warning: Unknown input " << line[0] << "!" << endl;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  //*********** Fill in the default values if not specified ***********//
  fillDefaults();

  //*********** Verify inputs ***********//
  verifyInputs();
  printf("%-40s %-s\n\n", "Finished Reading Input File:", fileName);
}

void ConfigSetup::fillDefaults(void)
{
  if(sys.elect.ewald == true) {
    sys.elect.enable = true;
  }

  if(sys.moves.intraSwap == DBL_MAX) {
    sys.moves.intraSwap = 0.000;
    printf("%-40s %-4.4f \n", "Default: Intra-Swap move frequency",
           sys.moves.intraSwap);
  }

  if(sys.moves.multiParticle == DBL_MAX) {
    sys.moves.multiParticle = 0.000;
    printf("%-40s %-4.4f \n",
           "Default: Multi-Particle move frequency",
           sys.moves.multiParticle);
  }

  if(sys.exclude.EXCLUDE_KIND == UINT_MAX) {
    sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONEFOUR_KIND;
    printf("%-40s %-s \n", "Default: Exclude", "ONE-FOUR");
  }

  if(sys.elect.enable && sys.elect.oneFourScale == DBL_MAX) {
    sys.elect.oneFourScale = 0.0f;
    if(sys.exclude.EXCLUDE_KIND != sys.exclude.EXC_ONEFOUR_KIND) {
      printf("%-40s %-lf \n", "Default: Modified 1-4 Electrostatic scaling",
             sys.elect.oneFourScale);
    }
  }

  if(in.prng.kind == "") {
    in.prng.kind = in.prng.KIND_RANDOM;
    printf("%-40s %-s \n", "Default: Random seed", "Active");
  }

#if ENSEMBLE == GEMC
  if(sys.gemc.kind == UINT_MAX) {
    sys.gemc.kind = mv::GEMC_NVT;
    printf("Default: Running NVT_GEMC\n");
  }
#endif

  if (sys.elect.ewald == true && sys.elect.readCache == false) {
    sys.elect.cache = true;
    printf("%-40s %-s \n", "Default: Cache Ewald Fourier", "Active");
  }

  if(sys.elect.enable && sys.elect.dielectric == DBL_MAX && in.ffKind.isMARTINI) {
    sys.elect.dielectric = 15.0f;
    printf("%-40s %-4.4f \n", "Default: Dielectric", sys.elect.dielectric);
  }

  if(sys.ff.cutoffLow == DBL_MAX) {
    sys.ff.cutoffLow = 0.00;
    printf("%-40s %-4.4f \n", "Default: Short Range Cutoff", sys.ff.cutoffLow);
  }

  out.state.files.psf.name = out.statistics.settings.uniqueStr.val +
                             "_merged.psf";
  for(int i = 0; i < BOX_TOTAL; i++) {
    if(i == 0)
      out.state.files.pdb.name[0] = out.statistics.settings.uniqueStr.val +
                                    "_BOX_0.pdb";
    else if(i == 1)
      out.state.files.pdb.name[1] = out.statistics.settings.uniqueStr.val +
                                    "_BOX_1.pdb";
  }
  out.state.files.seed.name = out.statistics.settings.uniqueStr.val + ".dat";
}

void ConfigSetup::verifyInputs(void)
{
  int i;
  if(!sys.elect.enable && sys.elect.oneFourScale != DBL_MAX) {
    printf("Warning: 1-4 Electrostatic scaling set, but will be ignored.\n");
  }

  if((sys.elect.oneFourScale != 0.0) &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONEFOUR_KIND)) {
    printf("Warning: 1-4 Electrostatic scaling set, but will be ignored.\n");
    sys.elect.oneFourScale = 0.0f;
  }

  if(sys.elect.ewald == false && sys.elect.enable == true) {
    printf("%-40s %-s \n",
           "Warning: Electrostatic calculation with Ewald method", "Inactive");
  }

  if(in.restart.enable  && sys.volume.hasVolume) {
    printf("Warning: Cell dimension set, but will be ignored in restart mode.\n");
  }

  if(in.prng.kind == "RANDOM" && in.prng.seed != UINT_MAX) {
    printf("Warning: Seed value set, but will be ignored.\n");
  }

  // Set output files
  if(out.statistics.settings.uniqueStr.val == "") {
    std::cout << "Error: Output name is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }

#if ENSEMBLE == GEMC
  if(sys.gemc.kind == mv::GEMC_NPT && sys.gemc.pressure == DBL_MAX) {
    std::cout << "Error: Pressure is not specified for NPT-GEMC!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.gemc.kind == mv::GEMC_NVT && sys.gemc.pressure != DBL_MAX) {
    std::cout << "Warning: Input pressure set, but will be ignored in NVT-GEMC"
              << std::endl;
  }
#endif
#if ENSEMBLE == NPT
  if(sys.gemc.pressure == DBL_MAX) {
    std::cout << "Error: Pressure is not specified for NPT!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.volume.cstVolBox0) {
    std::cout << "Warning: Fix volume of box 1 set, but will be ignored.\n";
    exit(EXIT_FAILURE);
  }
#endif

  if(in.prng.kind == "INTSEED" && in.prng.seed == UINT_MAX) {
    std::cout << "Error: Seed value is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(in.ffKind.numOfKinds == 0) {
    std::cout << "Error: Force field type is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(in.ffKind.numOfKinds > 1) {
    std::cout << "Error: Multiple Parameter types are specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if((!in.ffKind.isMARTINI && !in.ffKind.isEXOTIC) &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETWO_KIND)) {
    std::cout << "Warning: Exclude 1-2 is set for CHARMM type parameter.\n";
  }
  if(in.ffKind.isEXOTIC &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETWO_KIND)) {
    std::cout << "Warning: Exclude 1-2 is set for EXOTIC type parameter.\n";
  }
  if(in.ffKind.isEXOTIC &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETHREE_KIND)) {
    std::cout << "Warning: Exclude 1-3 is set for EXOTIC type parameter.\n";
  }
  if(in.files.param.name == "") {
    std::cout << "Error: Parameter file name is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.ff.VDW_KIND == UINT_MAX) {
    std::cout << "Error: Potential type is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND && sys.ff.doTailCorr == false) {
    std::cout << "Warning: Long Range Correction is Inactive for " <<
              "Non-truncated potential." << std::endl;
  }
  if(((sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND) ||
      (sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND)) && sys.ff.doTailCorr) {
    std::cout << "Warning: Long Range Correction is Active for " <<
              "truncated potential." << std::endl;
  }
  if(sys.ff.cutoff == DBL_MAX) {
    std::cout << "Error: Cutoff is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(sys.elect.ewald && (sys.elect.tolerance == DBL_MAX)) {
    std::cout << "Error: Tolerance is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.step.adjustment == ULONG_MAX) {
    std::cout << "Error: Move adjustment frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.step.equil == ULONG_MAX) {
    std::cout << "Error: Equilibration steps is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.step.total == ULONG_MAX) {
    std::cout << "Error: Total run steps is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.step.adjustment > sys.step.equil) {
    std::cout << "Error: Move adjustment frequency should be smaller " <<
              "than Equilibration steps!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.step.equil > sys.step.total) {
    std::cout << "Error: Equilibration steps should be smaller than " <<
              "Total run steps!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.moves.displace == DBL_MAX) {
    std::cout << "Error: Displacement move frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.moves.multiParticle == DBL_MAX) {
    std::cout <<
      "Error: Multi-Particle move frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.moves.rotate == DBL_MAX) {
    std::cout << "Error: Rotation move frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.moves.intraSwap == DBL_MAX) {
    std::cout << "Error: Intra-Swap move frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
#if ENSEMBLE == NPT
  if(sys.moves.volume == DBL_MAX) {
    std::cout << "Error: Volume move frequency is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif
#if ENSEMBLE == GEMC
  if(sys.moves.volume == DBL_MAX) {
    std::cout << "Error: Volume move frequency is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.moves.transfer == DBL_MAX) {
    std::cout << "Error: Molecule swap move frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(abs(sys.moves.displace + sys.moves.multiParticle +
         sys.moves.rotate + sys.moves.transfer +
         sys.moves.intraSwap + sys.moves.volume - 1.0) > 0.01) {
    std::cout << "Error: Sum of move frequncies are not equal to one!\n";
    exit(EXIT_FAILURE);
  }
#elif ENSEMBLE == NPT
  if(sys.moves.volume == DBL_MAX) {
    std::cout << "Error: Volume move frequency is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(abs(sys.moves.displace + sys.moves.multiParticle +
         sys.moves.rotate + sys.moves.intraSwap +
         sys.moves.volume - 1.0) > 0.01) {
    std::cout << "Error: Sum of move frequncies are not equal to one!\n";
    exit(EXIT_FAILURE);
  }

#elif ENSEMBLE == GCMC
  if(sys.moves.transfer == DBL_MAX) {
    std::cout << "Error: Molecule swap move frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(abs(sys.moves.displace + sys.moves.multiParticle + 
         sys.moves.rotate + sys.moves.intraSwap +
         sys.moves.transfer - 1.0) > 0.01) {
    std::cout << "Error: Sum of move frequncies are not equal to one!!\n";
    exit(EXIT_FAILURE);
  }
#else
  if(abs(sys.moves.displace + sys.moves.multiParticle + 
         sys.moves.rotate + sys.moves.intraSwap - 1.0) > 0.01) {
    std::cout << "Error: Sum of move frequncies are not equal to one!!\n";
    exit(EXIT_FAILURE);
  }
#endif

  for(i = 0 ; i < BOX_TOTAL ; i++) {
    if(in.files.pdb.name[i] == "") {
      std::cout << "Error: PDB file is not been specified for box number "
                << i << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  for(i = 0 ; i < BOX_TOTAL ; i++) {
    if(in.files.psf.name[i] == "") {
      std::cout << "Error: PSF file is not specified for box number " <<
                i << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if(!sys.volume.hasVolume && !in.restart.enable) {
    std::cout << "Error: This simulation requires to define " << 3 * BOX_TOTAL <<
              " cell basis vectors!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND && sys.ff.rswitch == DBL_MAX) {
    std::cout << "Error: Switch distance is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(((sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND) ||
      (sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND)) && sys.ff.rswitch != DBL_MAX) {
    std::cout << "Warning: Switch distance set, but will be ignored."
              << std::endl;
  }
  if(sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND &&
      sys.ff.rswitch >= sys.ff.cutoff) {
    std::cout << "Error: Switch distance should be less than Cutoff!\n";
    exit(EXIT_FAILURE);
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if(sys.cbmcTrials.bonded.ang == UINT_MAX) {
    std::cout << "Error: CBMC number of angle trials is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.cbmcTrials.bonded.dih == UINT_MAX) {
    std::cout << "Error: CBMC number of dihedral trials is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.cbmcTrials.nonbonded.first == UINT_MAX) {
    std::cout << "Error: CBMC number of first site trials is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(sys.cbmcTrials.nonbonded.nth == UINT_MAX) {
    std::cout << "Error: CBMC number of nth site trials is not specified!\n";
    exit(EXIT_FAILURE);
  }
#endif
  if(sys.T.inKelvin == DBL_MAX) {
    std::cout << "Error: Temperature is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.statistics.settings.uniqueStr.val == "") {
    std::cout << "Error: Output name is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.console.enable && out.console.frequency == ULONG_MAX) {
    std::cout << "Error: Console output frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.restart.settings.enable && out.restart.settings.frequency == ULONG_MAX) {
    std::cout << "Error: Restart coordinate frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.state.settings.enable && out.state.settings.frequency == ULONG_MAX) {
    std::cout << "Error: Coordinate frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.statistics.settings.block.enable &&
      out.statistics.settings.block.frequency == ULONG_MAX) {
    std::cout << "Error: Average output frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
#if ENSEMBLE == GCMC
  if(out.statistics.settings.hist.enable &&
      out.statistics.settings.hist.frequency == ULONG_MAX) {
    std::cout << "Error: Histogram output frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.state.files.hist.histName == "") {
    std::cout << "Error: Distribution file name is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(out.state.files.hist.sampleName == "") {
    std::cout << "Error: Histogram file name of is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(out.state.files.hist.letter == "") {
    std::cout << "Error: Run Letter of histogram file name is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.state.files.hist.number == "") {
    std::cout << "Error: Run number of histogram file is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if(out.state.files.hist.stepsPerHistSample == UINT_MAX) {
    std::cout << "Error: Histogram output sample frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
#endif
  if(!out.statistics.settings.block.enable && out.statistics.vars.energy.block) {
    printf("Note: Average output Inactived. Energy average output will be ignored.\n");
    out.statistics.vars.energy.block = false;
  }
  if(!out.statistics.settings.block.enable &&
      out.statistics.vars.pressure.block) {
    printf("Note: Average output Inactived. Pressure average output will be ignored.\n");
    out.statistics.vars.pressure.block = false;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.pressure.block) {
    printf("Note: Pressure Calculation Inactived. Pressure average output will be ignored.\n");
    out.statistics.vars.pressure.block = false;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.surfaceTension.block) {
    printf("Note: Pressure Calculation Inactived. Surface Tension average output will be ignored.\n");
    out.statistics.vars.surfaceTension.block = false;
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if(!out.statistics.settings.block.enable && out.statistics.vars.molNum.block) {
    printf("Note: Average output Inactived. Molecule average output will be ignored.\n");
    out.statistics.vars.molNum.block = false;
  }
#endif
  if(!out.statistics.settings.block.enable && out.statistics.vars.density.block) {
    printf("Note: Average output Inactived. Density average output will be ignored.\n");
    out.statistics.vars.density.block = false;
  }
#ifdef VARIABLE_VOLUME
  if(!out.statistics.settings.block.enable && out.statistics.vars.volume.block) {
    printf("Note: Average output Inactived. Volume average output will be ignored.\n");
    out.statistics.vars.volume.block = false;
  }
#endif
  if(!out.console.enable && out.statistics.vars.energy.fluct) {
    printf("Note: Console output Inactived. Energy output will be ignored.\n");
    out.statistics.vars.energy.fluct = false;
  }
  if(!out.console.enable && out.statistics.vars.pressure.fluct) {
    printf("Note: Console output Inactived. Pressure output will be ignored.\n");
    out.statistics.vars.pressure.fluct = false;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.pressure.fluct) {
    printf("Note: Pressure Calculation Inactived. Pressure output will be ignored.\n");
    out.statistics.vars.pressure.fluct = false;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.surfaceTension.fluct) {
    printf("Note: Pressure Calculation Inactived. Surface Tension output will be ignored.\n");
    out.statistics.vars.surfaceTension.fluct = false;
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if(!out.console.enable && out.statistics.vars.molNum.fluct) {
    printf("Note: Console output Inactived. Molecule output will be ignored.\n");
  }
#endif
  if(!out.console.enable && out.statistics.vars.density.fluct) {
    printf("Note: Console output Inactived. Density output will be ignored.\n");
  }
#ifdef VARIABLE_VOLUME
  if(!out.console.enable && out.statistics.vars.volume.fluct) {
    printf("Note: Console output Inactived. Volume output will be ignored.\n");
  }
#endif
}

const std::string config_setup::PRNGKind::KIND_RANDOM = "RANDOM",
                                          config_setup::PRNGKind::KIND_SEED = "INTSEED",
                                                                  config_setup::PRNGKind::KIND_RESTART = "RESTART",
                                                                                          config_setup::FFKind::FF_CHARMM = "CHARMM",
                                                                                                                config_setup::FFKind::FF_EXOTIC = "EXOTIC",
                                                                                                                                      config_setup::FFKind::FF_MARTINI = "MARTINI",
                                                                                                                                                            config_setup::FFValues::VDW = "VDW",
                                                                                                                                                                                    config_setup::FFValues::VDW_SHIFT = "VDW_SHIFT",
                                                                                                                                                                                                            config_setup::FFValues::VDW_SWITCH = "VDW_SWITCH",
                                                                                                                                                                                                                                    config_setup::Exclude::EXC_ONETWO = "1-2",
                                                                                                                                                                                                                                                           config_setup::Exclude::EXC_ONETHREE = "1-3",
                                                                                                                                                                                                                                                                                  config_setup::Exclude::EXC_ONEFOUR = "1-4";

const char ConfigSetup::defaultConfigFileName[] = "in.dat";
const char ConfigSetup::configFileAlias[] = "GO-MC Configuration File";

const uint config_setup::FFValues::VDW_STD_KIND = 0,
                                   config_setup::FFValues::VDW_SHIFT_KIND = 1,
                                                           config_setup::FFValues::VDW_SWITCH_KIND = 2,
                                                                                   config_setup::Exclude::EXC_ONETWO_KIND = 0,
                                                                                                          config_setup::Exclude::EXC_ONETHREE_KIND = 1,
                                                                                                                                 config_setup::Exclude::EXC_ONEFOUR_KIND = 2;
