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
  for(i=0; i<BOX_TOTAL; i++)
  {
    in.files.pdb.name[i] == "";
    in.files.psf.name[i] == "";
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
  out.statistics.settings.uniqueStr.val == "";
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
  for(k = 0; k < str.length(); k++)
  {
    str[k] = toupper(str[k]);
  }

  if(str == "ON" || str == "TRUE" || str == "YES")
    return true;
  else if(str == "OFF" || str == "FALSE" || str == "NO")
    return false;
  std::cout << "Error: " << str << "couldn't be recognized!" << std::endl;
  exit(0);
}

void ConfigSetup::Init(const char *fileName)
{
  std::vector<std::string> line;

  reader.Open(fileName);
  while(reader.readNextLine(line))
  {
    if(line[0] == "Restart")
    {
      in.restart.enable = checkBool(line[1]);
      printf("%-30s %-s \n", "Info: Restart simulation",  "Active");
    }
    else if(line[0] == "FirstStep")
    {
      in.restart.step = stringtoi(line[1]);
    }
    else if(line[0] == "PRNG")
    {
      in.prng.kind = line[1];
      printf("Info: Random seed Active.\n");
    }
    else if(line[0] == "Random_Seed")
    {
      in.prng.seed = stringtoi(line[1]);
      printf("Info: Constant seed Active.\n");
    }
    else if(line[0] == "ParaTypeCHARMM")
    {
      if(checkBool(line[1]))
      {
        in.ffKind.numOfKinds++;
        in.ffKind.isEXOTIC = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isCHARMM = true;
	printf("Info: PARAMETER file: CHARMM format!\n");
      }
    }
    else if(line[0] == "ParaTypeEXOTIC")
    {
      if(checkBool(line[1]))
      {
        in.ffKind.numOfKinds++;
        in.ffKind.isCHARMM = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isEXOTIC = true;
	printf("Info: PARAMETER file: EXOTIC format!\n");
      }
    }
    else if(line[0] == "ParaTypeMARTINI")
    {
      if(checkBool(line[1]))
      {
        in.ffKind.numOfKinds ++;
        in.ffKind.isEXOTIC = false;
        in.ffKind.isMARTINI = true;
        in.ffKind.isCHARMM = true;
	printf("Info: PARAMETER file: MARTINI using CHARMM format!\n");
      }
    }
    else if(line[0] == "Parameters")
    {
      in.files.param.name = line[1];
    }
    else if(line[0] == "Coordinates")
    {
      uint boxnum = stringtoi(line[1]);
      if(boxnum >= BOX_TOTAL)
      {
        std::cout<< "Error: Simulation requires " << BOX_TOTAL << " PDB file(s)!" << std::endl;
        exit(0);
      }
      in.files.pdb.name[boxnum] = line[2];
    }
    else if(line[0] == "Structure")
    {
      uint boxnum = stringtoi(line[1]);
      if(boxnum >= BOX_TOTAL)
      {
        std::cout<< "Error: Simulation requires " << BOX_TOTAL << " PSF file(s)!" << std::endl;
        exit(0);
      }
      in.files.psf.name[boxnum] = line[2];
    }
#if ENSEMBLE == GEMC
    else if(line[0] == "GEMC")
    {
      if(line[1] == "NVT")
      {
        sys.gemc.kind = mv::GEMC_NVT;
	printf("Info: Running NVT_GEMC.\n");
      }
      else if(line[1] == "NPT")
      {
        sys.gemc.kind = mv::GEMC_NPT;
	printf("Info: Running NPT_GEMC.\n");
      }
    }
    else if(line[0] == "Pressure")
    {
      sys.gemc.pressure = stringtod(line[1]);
      printf("%-30s %-4.4f bar.\n", "Info: Input Pressure", sys.gemc.pressure);
      sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
    }
#endif
#if ENSEMBLE == NPT
    else if(line[0] == "Pressure")
    {
      sys.gemc.kind = mv::GEMC_NPT;
      sys.gemc.pressure = stringtod(line[1]);
      printf("%-30s %-4.4f bar.\n", "Info: Input Pressure", sys.gemc.pressure);
      sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
    }
#endif
    else if(line[0] == "Temperature")
    {
      sys.T.inKelvin = stringtod(line[1]);
      printf("%-30s %-4.4f K.\n", "Info: Input Temperature", sys.T.inKelvin);
    }
    else if(line[0] == "Potential")
    {
      if(line[1] == "VDW")
      {
        sys.ff.VDW_KIND = sys.ff.VDW_STD_KIND;
	printf("%-30s %-s \n", "Info: Non-truncated potential", "Active");
      }
      else if(line[1] == "SHIFT")
      {
        sys.ff.VDW_KIND = sys.ff.VDW_SHIFT_KIND;
        printf("%-30s %-s \n", "Info: Shift truncated potential", "Active");
      }
      else if(line[1] == "SWITCH")
      {
        sys.ff.VDW_KIND = sys.ff.VDW_SWITCH_KIND;
	printf("%-30s %-s \n", "Info: Switch truncated potential", "Active");
      }
    }
    else if(line[0] == "LRC")
    {
      sys.ff.doTailCorr = checkBool(line[1]);
    }
    else if(line[0] == "Rswitch")
    {
      sys.ff.rswitch = stringtod(line[1]);
    }
    else if(line[0] == "Rcut")
    {
      sys.ff.cutoff = stringtod(line[1]);
      printf("%-30s %-4.4f A.\n", "Info: Cutoff", sys.ff.cutoff);
    }
    else if(line[0] == "RcutLow")
    {
      sys.ff.cutoffLow = stringtod(line[1]);
      printf("%-30s %-4.4f A.\n", "Info: Short Range Cutoff", sys.ff.cutoffLow);
    }
    else if(line[0] == "Exclude")
    {
      if(line[1] == sys.exclude.EXC_ONETWO)
      {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETWO_KIND;
	printf("%-30s %-s \n", "Info: Exclude", "ONE-TWO");
      }
      else if(line[1] == sys.exclude.EXC_ONETHREE)
      {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETHREE_KIND;
	printf("%-30s %-s \n", "Info: Exclude", "ONE-THREE");
      }
      else if(line[1] == sys.exclude.EXC_ONEFOUR)
      {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONEFOUR_KIND;
	printf("%-30s %-s \n", "Info: Exclude", "ONE-FOUR");
      }
    }
    else if(line[0] == "Ewald")
    {
      sys.elect.ewald = checkBool(line[1]);
      sys.elect.readEwald = true;
      if(sys.elect.ewald)
      {
	printf("%-30s %-s \n", "Info: Ewald Summation" , "Active");
      }
    }
    else if(line[0] == "ElectroStatic")
    {
      sys.elect.enable = checkBool(line[1]);
      sys.elect.readElect = true;
    }
    else if(line[0] == "Tolerance")
    {
      sys.elect.tolerance = stringtod(line[1]);
      sys.elect.alpha = sqrt(-1 * log(sys.elect.tolerance))/
                        sys.ff.cutoff;
      sys.elect.recip_rcut = 2 * (-log(sys.elect.tolerance))/
                             sys.ff.cutoff;
      printf("%-30s %-1.3E \n", "Info: Ewald Summation Tolerance" ,
	     sys.elect.tolerance);
    }
    else if(line[0] == "CachedFourier")
    {
      sys.elect.cache = checkBool(line[1]);
      sys.elect.readCache = true;
      if(sys.elect.cache)
      {
	printf("%-30s %-s \n", "Info: Cache Ewald Fourier", "Active");
      }
      else
      {
	printf("%-30s %-s \n", "Info: Cache Ewald Fourier", "Deactive");
      }
    }
    else if(line[0] == "1-4scaling")
    {
      sys.elect.oneFourScale = stringtod(line[1]);
    }
    else if(line[0] == "Dielectric")
    {
      sys.elect.dielectric = stringtod(line[1]);
      printf("%-30s %-4.4f \n", "Info: Dielectric", sys.elect.dielectric);
    }
    else if(line[0] == "RunSteps")
    {
      sys.step.total = stringtoi(line[1]);
      printf("%-30s %-d \n", "Info: Total number of steps", sys.step.total);
    }
    else if(line[0] == "EqSteps")
    {
      sys.step.equil = stringtoi(line[1]);
      printf("%-30s %-d \n", "Info: Number of equilibration steps",
	     sys.step.equil);
    }
    else if(line[0] == "AdjSteps")
    {
      sys.step.adjustment = stringtoi(line[1]);
      printf("%-30s %-d \n", "Info: Move adjustment frequency",
	     sys.step.adjustment);
    }
    else if(line[0] == "PressureCalc")
    {
      sys.step.pressureCalc = checkBool(line[1]);

      if(sys.step.pressureCalc && (line.size() == 3))
      {
	sys.step.pressureCalcFreq = stringtoi(line[2]);
	printf("%-30s %-d \n", "Info: Pressure calculation frequency",
	       sys.step.pressureCalcFreq);
      }
      else if(sys.step.pressureCalc && (line.size() == 2))
      {
	std::cout << "Error: Pressure calculation frequency has not been set!\n";
	exit(0);
      }
      else if(!sys.step.pressureCalc)
      {
        printf("%-30s %-s \n", "Info: Pressure calculation", "Deactive");
      }
    }
    else if(line[0] == "DisFreq")
    {
      sys.moves.displace = stringtod(line[1]);
      printf("%-30s %-4.4f \n", "Info: Displacement move frequency",
	     sys.moves.displace);
    }
    else if(line[0] == "IntraSwapFreq")
    {
      sys.moves.intraSwap = stringtod(line[1]);
      printf("%-30s %-4.4f \n", "Info: Intra-Swap move frequency",
	     sys.moves.intraSwap);
    }
    else if(line[0] == "RotFreq")
    {
      sys.moves.rotate = stringtod(line[1]);
      printf("%-30s %-4.4f \n", "Info: Rotation move frequency",
	     sys.moves.rotate);
    }
#ifdef VARIABLE_VOLUME
    else if(line[0] == "VolFreq")
    {
      sys.moves.volume = stringtod(line[1]);
      printf("%-30s %-4.4f \n", "Info: Volume move frequency",
	     sys.moves.volume);
    }
    else if(line[0] == "useConstantArea")
    {
      sys.volume.cstArea = checkBool(line[1]);
      if(sys.volume.cstArea)
	printf("Info: Volume change using constant X-Y area.\n");
      else
	printf("Info: Volume change using constant ratio.\n");
    }
    else if(line[0] == "FixVolBox0")
    {
      sys.volume.cstVolBox0 = checkBool(line[1]);
      if (sys.volume.cstVolBox0)
        printf("%-30s %-d \n", "Info: Fix volume box", 1);
    }
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
    else if(line[0] == "SwapFreq")
    {
#if ENSEMBLE == NVT || ENSEMBLE == NPT
      sys.moves.transfer = 0.000;
#else
      sys.moves.transfer = stringtod(line[1]);
#endif
      printf("%-30s %-4.4f \n", "Info: Molecule swap  move frequency",
	     sys.moves.transfer);
    }
#endif
    else if(line[0] == "BoxDim")
    {
      uint box = stringtoi(line[1]);
      if(box < BOX_TOTAL)
      {
        XYZ temp;
        sys.volume.boxCoordRead++;
        sys.volume.hasVolume = (sys.volume.boxCoordRead == BOX_TOTAL);
        temp.x = stringtod(line[2]);
        temp.y = stringtod(line[3]);
        temp.z = stringtod(line[4]);
        sys.volume.axis.Set(box, temp);
	printf("%-30s: %-5d %-4.4f %-4.4f %-4.4f \n", "Info: Simulation dimension of box", box, temp.x, temp.y, temp.z);
      }
      else
      {
        std::cout<< "Error: This simulation requires only " << BOX_TOTAL << " number of box dimension(s)!" << std::endl;
        exit(0);
      }
    }
#ifdef VARIABLE_PARTICLE_NUMBER
    else if(line[0] == "CBMC_First")
    {
      sys.cbmcTrials.nonbonded.first = stringtoi(line[1]);
    }
    else if(line[0] == "CBMC_Nth")
    {
      sys.cbmcTrials.nonbonded.nth = stringtoi(line[1]);
    }
    else if(line[0] == "CBMC_Ang")
    {
      sys.cbmcTrials.bonded.ang = stringtoi(line[1]);
    }
    else if(line[0] == "CBMC_Dih")
    {
      sys.cbmcTrials.bonded.dih = stringtoi(line[1]);
    }
#endif
#if ENSEMBLE == GCMC
    else if(line[0] == "ChemPot")
    {
      if(line.size() != 3)
      {
        std::cout << "Error: Chemical potential parameters have not been specified!" << std::endl;
        exit(0);
      }
      std::string resName = line[1];
      double val = stringtod(line[2]);
      sys.chemPot.cp[resName] = val;
      printf("%-30s %-6s %-6.4f K.\n", "Info: Chemical potential:",
	     resName.c_str(), val);
    }
    else if(line[0] == "Fugacity")
    {
      if(line.size() != 3)
      {
	std::cout << "Error: Fugacity parameters have not been specified!" << std::endl;
	exit(0);
      }
      sys.chemPot.isFugacity = true;
      std::string resName = line[1];
      double val = stringtod(line[2]);
      sys.chemPot.cp[resName] = val * unit::BAR_TO_K_MOLECULE_PER_A3;
      printf("%-30s %-6s %-6.4f bar.\n", "Info: Fugacity:", resName.c_str(),
	     val);
    }
#endif
    else if(line[0] == "OutputName")
    {
      out.statistics.settings.uniqueStr.val = line[1];
      printf("%-30s %-s \n", "Info: Output name", line[1].c_str());
    }
    else if(line[0] == "CoordinatesFreq")
    {
      out.state.settings.enable = checkBool(line[1]);
      out.state.settings.frequency = stringtoi(line[2]);
      if(out.state.settings.enable)
      {
	printf("%-30s %-d \n", "Info: Coordinate frequency:",
	       out.state.settings.frequency);
      }
      else
	printf("%-30s %-s \n", "Info: Printing coordinate", "Deactive");
    }
    else if(line[0] == "RestartFreq")
    {
      out.restart.settings.enable = checkBool(line[1]);
      out.restart.settings.frequency = stringtoi(line[2]);
      if(out.restart.settings.enable)
      {
	printf("%-30s %-d \n", "Info: Restart frequency:",
	       out.restart.settings.frequency);
      }
      else
	printf("%-30s %-s \n", "Info: Printing restart coordinate", "Deactive");
    }
    else if(line[0] == "ConsoleFreq")
    {
      out.console.enable = checkBool(line[1]);
      out.console.frequency = stringtoi(line[2]);
      if(out.console.enable)
      {
	printf("%-30s %-d \n", "Info: Console output frequency:",
	       out.console.frequency);
      }
      else
	printf("%-30s %-s \n", "Info: Console output", "Deactive");
    }
    else if(line[0] == "BlockAverageFreq")
    {
      out.statistics.settings.block.enable = checkBool(line[1]);
      out.statistics.settings.block.frequency = stringtoi(line[2]);
      if(out.statistics.settings.block.enable)
      {
	printf("%-30s %-d \n", "Info: Average output frequency:",
	       out.statistics.settings.block.frequency);
      }
      else
	printf("%-30s %-s \n", "Info: Average output", "Deactive");
    }
#if ENSEMBLE == GCMC
    else if(line[0] == "HistogramFreq")
    {
      out.statistics.settings.hist.enable = checkBool(line[1]);
      out.statistics.settings.hist.frequency = stringtoi(line[2]);
      if(out.statistics.settings.hist.enable)
      {
	printf("%-30s %-d \n", "Info: Histogram output frequency:",
	       out.statistics.settings.hist.frequency);
      }
      else
	printf("%-30s %-s \n", "Info: Histogram output", "Deactive");
    }
    else if(line[0] == "DistName")
    {
      out.state.files.hist.histName = line[1];
    }
    else if(line[0] == "HistName")
    {
      out.state.files.hist.sampleName = line[1];
    }
    else if(line[0] == "RunNumber")
    {
      out.state.files.hist.number = line[1];
    }
    else if(line[0] == "RunLetter")
    {
      out.state.files.hist.letter = line[1];
    }
    else if(line[0] == "SampleFreq")
    {
      out.state.files.hist.stepsPerHistSample = stringtoi(line[1]);
    }
#endif
    else if(line[0] == "OutEnergy")
    {
      out.statistics.vars.energy.block = checkBool(line[1]);
      out.statistics.vars.energy.fluct = checkBool(line[2]);
    }
    else if(line[0] == "OutPressure")
    {
      out.statistics.vars.pressure.block = checkBool(line[1]);
      out.statistics.vars.pressure.fluct = checkBool(line[2]);
    }
#ifdef VARIABLE_PARTICLE_NUMBER
    else if(line[0] == "OutMolNum")
    {
      out.statistics.vars.molNum.block = checkBool(line[1]);
      out.statistics.vars.molNum.fluct = checkBool(line[2]);
    }
#endif
    else if(line[0] == "OutDensity")
    {
      out.statistics.vars.density.block = checkBool(line[1]);
      out.statistics.vars.density.fluct = checkBool(line[2]);
    }
    else if(line[0] == "OutSurfaceTension")
    {
      out.statistics.vars.surfaceTension.block = checkBool(line[1]);
      out.statistics.vars.surfaceTension.fluct = checkBool(line[2]);
    }
#ifdef VARIABLE_VOLUME
    else if(line[0] == "OutVolume")
    {
      out.statistics.vars.volume.block = checkBool(line[1]);
      out.statistics.vars.volume.fluct = checkBool(line[2]);
    }
#endif
    else
    {
      cout << "Warning: Unknown input " << line[0] << "!" << endl;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  //*********** Fill in the default values if not specified ***********//
  fillDefaults();

  //*********** Verify inputs ***********//
  verifyInputs();
}

void ConfigSetup::fillDefaults(void)
{
  if(sys.moves.intraSwap == DBL_MAX)
  {
    std::cout << "By default intra box swap frequency has been set to zero" << std::endl;
    sys.moves.intraSwap = 0.000;
  }

  if(sys.exclude.EXCLUDE_KIND == UINT_MAX)
  {
    std::cout << "Warning: By default value (1-3) for exclude is selected!" << std::endl;
    sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETHREE_KIND;
  }

  if(in.prng.kind == "")
  {
    std::cout << "Warning: By default, random seed has been selected!" << std::endl;
    in.prng.kind = in.prng.KIND_RANDOM;
  }

#if ENSEMBLE == GEMC
  if(sys.gemc.kind == UINT_MAX)
  {
    std::cout << "Warning: By default, GEMC-NVT has been selected!" << std::endl;
    sys.gemc.kind = mv::GEMC_NVT;
  }
#endif
  if(sys.elect.enable && sys.elect.oneFourScale == DBL_MAX)
  {
    std::cout << "Warning: 1-4 electro static scaling has been set to zero!" << std::endl;
    sys.elect.oneFourScale = 0.0f;
  }

  if (sys.elect.ewald == true)
  {
    sys.elect.enable = true;
  }

  if (sys.elect.ewald == true && sys.elect.readCache == false)
  {
    sys.elect.cache = true;
    std::cout << "Warning: By default, Fourier terms of ewald method will be cached!" << std::endl;
  }

  if (sys.elect.ewald == false && sys.elect.enable == true)
  {
    std::cout << "Warning: Ewald method would not be used to calculate electrostatic energy!" << std::endl;
  }

  if (sys.elect.ewald == false && sys.elect.enable == false)
  {
    std::cout << "Warning: Electrostatic energy would not be calculated!" << std::endl;
  }

  if(sys.elect.enable && sys.elect.dielectric == DBL_MAX && in.ffKind.isMARTINI)
  {
    std::cout << "Warning: Dielectric will be set to 15.0 for Martini forcefield!" << std::endl;
    sys.elect.dielectric = 15.0f;
  }

  // Set output files
  if(out.statistics.settings.uniqueStr.val == "")
  {
    std::cout << "Error: Output name is required!" << std::endl;
    exit(0);
  }
  out.state.files.psf.name = out.statistics.settings.uniqueStr.val + "_merged.psf";
  for(int i = 0; i<BOX_TOTAL; i++)
  {
    if(i==0)
      out.state.files.pdb.name[0] = out.statistics.settings.uniqueStr.val + "_BOX_0.pdb";
    else if(i==1)
      out.state.files.pdb.name[1] = out.statistics.settings.uniqueStr.val + "_BOX_1.pdb";
  }
  out.state.files.seed.name = out.statistics.settings.uniqueStr.val + ".dat";
}

void ConfigSetup::verifyInputs(void)
{
  int i;
#if ENSEMBLE == GEMC
  if(sys.gemc.kind == mv::GEMC_NPT && sys.gemc.pressure == DBL_MAX)
  {
    std::cout << "Error: Pressure has not been specified for GEMC-NPT!" << std::endl;
    exit(0);
  }
  if(sys.gemc.kind == mv::GEMC_NVT && sys.gemc.pressure != DBL_MAX)
  {
    std::cout << "Warning: Pressure will be ignored for GEMC-NVT!" << std::endl;
  }
#endif
#if ENSEMBLE == NPT
  if(sys.gemc.pressure == DBL_MAX)
  {
    std::cout << "Error: Pressure has not been specified for NPT simulation!" << std::endl;
    exit(0);
  }
  if(sys.volume.cstVolBox0)
  {
    std::cout << "Note: Volume cannot be fix for NPT simulation.\n";
    exit(0);
  }
#endif

  if(in.restart.enable == true && in.restart.step == ULONG_MAX)
  {
    //std::cout << "Error: Restart step is needed!" << std::endl;
    //exit(0);
  }
  if(in.restart.enable == false && in.restart.step != ULONG_MAX)
  {
    //std::cout << "Warning: Restart step will not be used!" << std::endl;
  }
  if(in.prng.kind == "RANDOM" && in.prng.seed != UINT_MAX)
  {
    std::cout << "Warning: Seed will not be used for RANDOM type seed!" << std::endl;
  }
  if(in.prng.kind == "INTSEED" && in.prng.seed == UINT_MAX)
  {
    std::cout << "Error: Seed is required for INTSEED type seed!" << std::endl;
    exit(0);
  }
  if(in.ffKind.numOfKinds == 0)
  {
    std::cout << "Error: Force field type has not been defined!" << std::endl;
    exit(0);
  }
  if(in.ffKind.numOfKinds > 1)
  {
    std::cout << "Error: One type of Parameter type should be set!" << std::endl;
    exit(0);
  }
  if((!in.ffKind.isMARTINI && !in.ffKind.isEXOTIC) && (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETWO_KIND))
  {
    std::cout << "Error: 1-3 interaction is not valid for CHARMM force field!" << std::endl;
    exit(0);
  }
  if(in.ffKind.isEXOTIC && (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETWO_KIND))
  {
    std::cout << "Error: 1-3 interaction is not valid for EXOTIC force field!" << std::endl;
    exit(0);
  }
  if(in.ffKind.isEXOTIC && (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETHREE_KIND))
  {
    std::cout << "Error: 1-4 interaction is not valid for EXOTIC force field!" << std::endl;
    exit(0);
  }
  if(in.files.param.name == "")
  {
    std::cout << "Error: Parameter file name has not been defined!" << std::endl;
    exit(0);
  }
  if(sys.ff.VDW_KIND == UINT_MAX)
  {
    std::cout << "Error: Potential type has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND && sys.ff.doTailCorr == false)
  {
    std::cout << "Warning: Long Range Correction has been disabled for standard VDW potential type!" << std::endl;
  }
  if((sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND || sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND) && sys.ff.doTailCorr)
  {
    std::cout << "Warning: Long Range Correction will be disabled for shift and switch potential!" << std::endl;
  }
  if(sys.ff.cutoff == DBL_MAX)
  {
    std::cout << "Error: Cut off is required!" << std::endl;
    exit(0);
  }
  if(sys.ff.cutoffLow == DBL_MAX)
  {
    std::cout << "Warning: Cut off Low is set to 1 A!" << std::endl;
    sys.ff.cutoffLow = 1.00;
  }
  if(sys.elect.ewald && (sys.elect.tolerance == DBL_MAX))
  {
    std::cout << "Error: Tolerance has not been specified for Ewald summation method!" << std::endl;
    exit(0);
  }
  if(sys.step.adjustment == ULONG_MAX)
  {
    std::cout << "Error: Adjustment steps has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.step.equil == ULONG_MAX)
  {
    std::cout << "Error: Equilibrium steps has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.step.total == ULONG_MAX)
  {
    std::cout << "Error: Total run steps has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.step.adjustment > sys.step.equil)
  {
    std::cout << "Error: Adjustment steps should be smaller than Equilibrium steps!" << std::endl;
    exit(0);
  }
  if(sys.step.equil > sys.step.total)
  {
    std::cout << "Error: Equilibrium steps should be smaller than Total run steps!" << std::endl;
    exit(0);
  }
  if(sys.moves.displace == DBL_MAX)
  {
    std::cout << "Error: Displacement frequency has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.moves.rotate == DBL_MAX)
  {
    std::cout << "Error: Rotation frequency has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.moves.intraSwap == DBL_MAX)
  {
    std::cout << "Error: IntraSwap frequency has not been specified!" << std::endl;
    exit(0);
  }
#if ENSEMBLE == NPT
  if(sys.moves.volume == DBL_MAX)
  {
    std::cout << "Error: Volume swap frequency has not been specified!" << std::endl;
    exit(0);
  }
#endif
#if ENSEMBLE == GEMC
  if(sys.moves.volume == DBL_MAX)
  {
    std::cout << "Error: Volume swap frequency has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.moves.transfer == DBL_MAX)
  {
    std::cout << "Error: Molecule swap frequency has not been specified!" << std::endl;
    exit(0);
  }
  if(abs(sys.moves.displace + sys.moves.rotate + sys.moves.transfer + sys.moves.intraSwap + sys.moves.volume - 1.0) > 0.01)
  {
    std::cout << "Error: Sum of move frequncies are not equal to one!" << std::endl;
    exit(0);
  }
#elif ENSEMBLE == NPT
  if(sys.moves.volume == DBL_MAX)
  {
    std::cout << "Error: Volume swap frequency has not been specified!" << std::endl;
    exit(0);
  }
  if(abs(sys.moves.displace + sys.moves.rotate + sys.moves.intraSwap + sys.moves.volume - 1.0) > 0.01)
  {
    std::cout << "Error: Sum of move frequncies are not equal to one!" << std::endl;
    exit(0);
  }

#elif ENSEMBLE == GCMC
  if(sys.moves.transfer == DBL_MAX)
  {
    std::cout << "Error: Molecule swap frequency has not been specified!" << std::endl;
    exit(0);
  }
  if(abs(sys.moves.displace + sys.moves.rotate + sys.moves.intraSwap + sys.moves.transfer - 1.0) > 0.01)
  {
    std::cout << "Error: Sum of move frequncies are not equal to one!" << std::endl;
    exit(0);
  }
#else
  if(abs(sys.moves.displace + sys.moves.rotate + sys.moves.intraSwap - 1.0) > 0.01)
  {
    std::cout << "Error: Sum of move frequncies are not equal to one!" << std::endl;
    exit(0);
  }
#endif

  for(i = 0 ; i < BOX_TOTAL ; i++)
  {
    if(in.files.pdb.name[i] == "")
    {
      std::cout << "Error: PDB file has not been defined for Box number " << i << "!" <<std::endl;
      exit(0);
    }
  }
  for(i = 0 ; i < BOX_TOTAL ; i++)
  {
    if(in.files.psf.name[i] == "")
    {
      std::cout << "Error: PSF file has not been defined for Box number " << i << "!" <<std::endl;
      exit(0);
    }
  }
  if(!sys.volume.hasVolume)
  {
    std::cout << "Error: This simulation type requires to define " << BOX_TOTAL << " box dimentions!" <<std::endl;
    exit(0);
  }
  if(sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND && sys.ff.rswitch == DBL_MAX)
  {
    std::cout << "Error: Switch cut off has not been defined!" << std::endl;
    exit(0);
  }
  if((sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND || sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND) && sys.ff.rswitch != DBL_MAX)
  {
    std::cout << "Warning: Switch cut off will be ignored!" << std::endl;
  }
  if(sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND && sys.ff.rswitch >= sys.ff.cutoff)
  {
    std::cout << "Error: Switch cut off should be smaller than cut off!" << std::endl;
    exit(0);
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if(sys.cbmcTrials.bonded.ang == UINT_MAX)
  {
    std::cout << "Error: CBMC number of angle trials has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.cbmcTrials.bonded.dih == UINT_MAX)
  {
    std::cout << "Error: CBMC number of dihedral trials has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.cbmcTrials.nonbonded.first == UINT_MAX)
  {
    std::cout << "Error: CBMC number of first site trials has not been specified!" << std::endl;
    exit(0);
  }
  if(sys.cbmcTrials.nonbonded.nth == UINT_MAX)
  {
    std::cout << "Error: CBMC number of nth site trials has not been specified!" << std::endl;
    exit(0);
  }
#endif
  if(sys.T.inKelvin == DBL_MAX)
  {
    std::cout << "Error: Temperature has not been specified!" << std::endl;
    exit(0);
  }
  if(out.statistics.settings.uniqueStr.val == "")
  {
    std::cout<< "Error: Output name has not been specified!" << std::endl;
    exit(0);
  }
  if(!out.state.settings.enable)
  {
    std::cout << "Warning: Coordinates output is disabled!" << std::endl;
  }
  if(!out.restart.settings.enable)
  {
    std::cout << "Warning: Restart coordinate output is disabled!" << std::endl;
  }
  if(!out.console.enable)
  {
    std::cout << "Warning: Console output frequency is disabled!" << std::endl;
  }
  if(!out.statistics.settings.block.enable)
  {
    std::cout << "Warning: Block average output is disabled!" << std::endl;
  }
  if(out.console.enable && out.console.frequency == ULONG_MAX)
  {
    if(sys.step.total > 1000)
    {
      out.console.frequency = (ulong)sys.step.total / 1000;
    }
    else
    {
      out.console.frequency = (ulong)sys.step.total / 100;
    }
    std::cout << "Warning: By default console output frequency has been set to " << out.console.frequency << "!" << std::endl;
  }
  if(out.restart.settings.enable && out.restart.settings.frequency == ULONG_MAX)
  {
    out.restart.settings.frequency = (ulong)sys.step.total;
    std::cout << "Warning: By default restart coordinate output frequency has been set to " << out.restart.settings.frequency << "!" << std::endl;
  }
  if(out.state.settings.enable && out.state.settings.frequency == ULONG_MAX)
  {
    out.state.settings.frequency = (ulong)sys.step.total / 10;
    std::cout << "Warning: By default coordinate output frequency has been set to " << out.state.settings.frequency << "!" << std::endl;
  }
  if(out.statistics.settings.block.enable && out.statistics.settings.block.frequency == ULONG_MAX)
  {
    out.statistics.settings.block.frequency = (ulong)sys.step.total / 100;
    std::cout << "Warning: By default block average output frequency has been set to " << out.statistics.settings.block.frequency << "!" << std::endl;
  }
  if(sys.step.pressureCalc && (sys.step.pressureCalcFreq == ULONG_MAX))
  {
    sys.step.pressureCalcFreq = (ulong)(out.statistics.settings.block.frequency / 100);
    std::cout << "Warning: By default pressure will be calculated every " << sys.step.pressureCalcFreq << " steps!";
  }
#if ENSEMBLE == GCMC
  if(!out.statistics.settings.hist.enable)
  {
    std::cout << "Warning: Histogram output is disabled!" << std::endl;
  }
  if(out.statistics.settings.hist.enable && out.statistics.settings.hist.frequency == ULONG_MAX)
  {
    if(sys.step.total > 1000)
    {
      out.statistics.settings.hist.frequency = (ulong)sys.step.total / 1000;
    }
    else
    {
      out.statistics.settings.hist.frequency = (ulong)sys.step.total / 100;
    }
    std::cout << "Warning: By default histogram output frequency has been set to " << out.statistics.settings.hist.frequency << "!" << std::endl;
  }
  if(out.state.files.hist.histName == "")
  {
    std::cout << "Error: Distribution file name has not been set!" << std::endl;
    exit(0);
  }
  if(out.state.files.hist.sampleName == "")
  {
    std::cout << "Error: Histogram file name of has not been set!" << std::endl;
    exit(0);
  }
  if(out.state.files.hist.letter == "")
  {
    std::cout << "Error: Run Letter of histogram file name has not been set!" << std::endl;
    exit(0);
  }
  if(out.state.files.hist.number == "")
  {
    std::cout << "Error: Run number of histogram file has not been set!" << std::endl;
    exit(0);
  }
  if(out.state.files.hist.stepsPerHistSample == UINT_MAX)
  {
    std::cout << "Error: Histogram output sample frequency has not been set!" << std::endl;
    exit(0);
  }
#endif
  if(!out.statistics.settings.block.enable && out.statistics.vars.energy.block)
  {
    std::cout<< "Warning: Output block average has been set off, energy output for block average will be ignored!" << std::endl;
  }
  if(!out.statistics.settings.block.enable && out.statistics.vars.pressure.block)
  {
    std::cout<< "Warning: Output block average has been set off, pressure output for block average will be ignored!" << std::endl;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.pressure.block)
  {
    std::cout<< "Warning: Pressure calculation has been turn off, pressure output for block average will be ignored!" << std::endl;
    out.statistics.vars.pressure.block = false;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.surfaceTension.block)
  {
    std::cout<< "Warning: Pressure calculation has been turn off, surface tension output for block average will be ignored!" << std::endl;
    out.statistics.vars.surfaceTension.block = false;
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if(!out.statistics.settings.block.enable && out.statistics.vars.molNum.block)
  {
    std::cout<< "Warning: Output block average has been set off, molecule number output for block average will be ignored!" << std::endl;
  }
#endif
  if(!out.statistics.settings.block.enable && out.statistics.vars.density.block)
  {
    std::cout<< "Warning: Output block average has been set off, density output for block average will be ignored!" << std::endl;
  }
#ifdef VARIABLE_VOLUME
  if(!out.statistics.settings.block.enable && out.statistics.vars.volume.block)
  {
    std::cout<< "Warning: Output block average has been set off, volume ouput for block average will be ignored!" << std::endl;
  }
#endif
  if(!out.console.enable && out.statistics.vars.energy.fluct)
  {
    std::cout<< "Warning: Console output has been set off, energy fluctuation ouput will be ignored!" << std::endl;
  }
  if(!out.console.enable && out.statistics.vars.pressure.fluct)
  {
    std::cout<< "Warning: Console output has been set off, pressure fluctuation ouput will be ignored!" << std::endl;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.pressure.fluct)
  {
    std::cout<< "Warning: Pressure calculation has been turn off, pressure fluctuation output will be ignored!" << std::endl;
    out.statistics.vars.pressure.fluct = false;
  }
  if(!sys.step.pressureCalc && out.statistics.vars.surfaceTension.fluct)
  {
    std::cout<< "Warning: Pressure calculation has been turn off, surface tension fluctuation output will be ignored!" << std::endl;
    out.statistics.vars.surfaceTension.fluct = false;
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if(!out.console.enable && out.statistics.vars.molNum.fluct)
  {
    std::cout<< "Warning: Console output has been set off, molecule number fluctuation ouput will be ignored!" << std::endl;
  }
#endif
  if(!out.console.enable && out.statistics.vars.density.fluct)
  {
    std::cout<< "Warning: Console output has been set off, density fluctuation ouput will be ignored!" << std::endl;
  }
#ifdef VARIABLE_VOLUME
  if(!out.console.enable && out.statistics.vars.volume.fluct)
  {
    std::cout<< "Warning: Console output has been set off, volume fluctuation ouput will be ignored!" << std::endl;
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
