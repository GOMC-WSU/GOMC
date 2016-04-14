/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.70 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <map> //for function handle storage.
#include <string> //for var names, etc.
#include <vector>
#include <string>

#include "ConfigSetup.h"
std::string prefix = ""; //DECLARED STRING PREFIX TO MAKE AVAILABLE GLOBALLY
//#define UINT_MAX 0xffffffff
//#define ULONG_MAX 0xffffffffUL
#define DBL_MAX 1.7976931348623158e+308

int stringtoi(const std::string& s) {
    std::istringstream str(s);
    int i;
    str >> i;
    return i;
}

double stringtod(const std::string& s) {
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
	sys.elect.tolerance = DBL_MAX;
	sys.elect.oneFourScale = DBL_MAX;
	sys.elect.dielectric = DBL_MAX;
	sys.step.total = ULONG_MAX;
	sys.step.equil = ULONG_MAX;
	sys.step.adjustment = ULONG_MAX;
	in.ffKind.numOfKinds = 0;
	sys.exclude.EXCLUDE_KIND = UINT_MAX;
	in.prng.kind = "";
	in.files.param.name = "";
	for(i=0;i<BOX_TOTAL;i++)
	{
		in.files.pdb.name[i] == "";
		in.files.psf.name[i] == "";
	}
#if ENSEMBLE == GEMC
	sys.gemc.kind = UINT_MAX;
	sys.gemc.pressure = DBL_MAX;
#endif
	sys.T.inKelvin = DBL_MAX;
	sys.ff.VDW_KIND = UINT_MAX;
	sys.ff.doTailCorr = true;
	sys.ff.rswitch = DBL_MAX;
	sys.ff.cutoff = DBL_MAX;
	sys.moves.displace = DBL_MAX;
	sys.moves.rotate = DBL_MAX;
	sys.moves.intraSwap = DBL_MAX;
	out.state.settings.enable = true;
	out.restart.settings.enable = true;
	out.console.enable = true;
	out.statistics.settings.fluct.enable = true;
	out.statistics.settings.block.enable = true;
#if ENSEMBLE == GCMC
	out.statistics.settings.hist.enable = true;
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
	out.statistics.settings.fluct.frequency = ULONG_MAX;
	out.statistics.vars.energy.block = false;
	out.statistics.vars.energy.fluct = false;
	out.statistics.vars.energy.hist = false;
	out.statistics.vars.pressure.block = false;
	out.statistics.vars.pressure.fluct = false;
	out.statistics.vars.pressure.hist = false;
#ifdef VARIABLE_PARTICLE_NUMBER
	sys.moves.transfer = DBL_MAX;
	sys.cbmcTrials.bonded.ang = UINT_MAX;
	sys.cbmcTrials.bonded.dih = UINT_MAX;
	sys.cbmcTrials.nonbonded.first = UINT_MAX;
	sys.cbmcTrials.nonbonded.nth = UINT_MAX;
	out.statistics.vars.molNum.block = false;
	out.statistics.vars.molNum.fluct = false;
	out.statistics.vars.molNum.hist = false;
	out.statistics.vars.acceptAngles.block = false;
	out.statistics.vars.acceptAngles.fluct = false;
	out.statistics.vars.acceptAngles.hist = false;
#endif
#ifdef VARIABLE_VOLUME
	sys.moves.volume = DBL_MAX;
	out.statistics.vars.volume.block = false;
	out.statistics.vars.volume.fluct = false;
	out.statistics.vars.volume.hist = false;
#endif
#ifdef VARIABLE_DENSITY
	out.statistics.vars.density.block = false;
	out.statistics.vars.density.fluct = false;
	out.statistics.vars.density.hist = false;
#endif
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

	//ADDED LINES TO GET PREFIX PATH AND 
	//CONVERT CHAR ARRY TO STRING TYPE TO GET PREFIX PATH
	//string prefix;
	std::string inputString(fileName);

	size_t found = inputString.rfind('/');
	if (found != std::string::npos)
		prefix = inputString.substr(0, inputString.find_last_of('/')) + '/';
	//END ADD LINE
	
	reader.Open(fileName);
	while(reader.readNextLine(line))
	{
		if(line[0] == "Restart") {
			in.restart.enable = checkBool(line[1]);
		} else if(line[0] == "FirstStep") {
			in.restart.step = stringtoi(line[1]);
		} else if(line[0] == "PRNG") {
			in.prng.kind = line[1];
		} else if(line[0] == "Random_Seed") {
			in.prng.seed = stringtoi(line[1]);
		} else if(line[0] == "ParaTypeCHARMM") {
			if(checkBool(line[1]))
			{
				in.ffKind.numOfKinds++;
				in.ffKind.isEXOTIC = false;
				in.ffKind.isMARTINI = false;
				in.ffKind.isCHARMM = true;
				std::cout << "REMINDER: CHARMM force field has been selected!" << std::endl;
			}
		} else if(line[0] == "ParaTypeEXOTIC") {
			if(checkBool(line[1]))
			{
				in.ffKind.numOfKinds++;
				in.ffKind.isCHARMM = false;
				in.ffKind.isMARTINI = false;
			        in.ffKind.isEXOTIC = true;
				std::cout << "REMINDER: EXOTIC force field has been selected!" << std::endl;
		}
		} else if(line[0] == "ParaTypeMARTINI") {
			if(checkBool(line[1]))
			{
				in.ffKind.numOfKinds ++;
				in.ffKind.isEXOTIC = false;
				in.ffKind.isMARTINI = true;
				in.ffKind.isCHARMM = true;
				std::cout << "REMINDER: MARTINI force field has been selected!" << std::endl;
			}
		} else if(line[0] == "Parameters") {
			in.files.param.name = prefix + line[1]; //ADDED PREFIX
		} else if(line[0] == "Coordinates") {
			uint boxnum = stringtoi(line[1]);
			if(boxnum >= BOX_TOTAL)
			{
				std::cout<< "Error: This simulation requires only " << BOX_TOTAL << " number of PDB file(s)!" << std::endl;
				exit(0);
			}
			in.files.pdb.name[boxnum] = prefix + line[2]; //ADDED PREFIX
		} else if(line[0] == "Structure") {
			uint boxnum = stringtoi(line[1]);
			if(boxnum >= BOX_TOTAL)
			{
				std::cout<< "Error: This simulation requires only " << BOX_TOTAL << " number of PSF file(s)!" << std::endl;
				exit(0);
			}
			in.files.psf.name[boxnum] = prefix + line[2]; //ADDED PREFIX
		}
#if ENSEMBLE == GEMC
		else if(line[0] == "GEMC")
		{
			if(line[1] == "NVT")
			{
				sys.gemc.kind = mv::GEMC_NVT;
				std::cout<< " NVT_GEMC simulation has been selected " << std::endl;
			}
			else if(line[1] == "NPT")
			{
				sys.gemc.kind = mv::GEMC_NPT;
				std::cout<< " NPT_GEMC simulation has been selected " << std::endl;
                        } 
		}
		else if(line[0] == "Pressure")
		{
			sys.gemc.pressure = stringtod(line[1]);
			std::cout<< " Pressure of system has been set to " << sys.gemc.pressure << " bar " << std::endl;
			sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
		}
#endif
		else if(line[0] == "Temperature")
		{
			sys.T.inKelvin = stringtod(line[1]);
		}
		else if(line[0] == "Potential")
		{
			if(line[1] == "VDW")
				sys.ff.VDW_KIND = sys.ff.VDW_STD_KIND;
			else if(line[1] == "SHIFT")
				sys.ff.VDW_KIND = sys.ff.VDW_SHIFT_KIND;
			else if(line[1] == "SWITCH")
				sys.ff.VDW_KIND = sys.ff.VDW_SWITCH_KIND;
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
		}
		else if(line[0] == "Exclude")
		{
			if(line[1] == sys.exclude.EXC_ONETWO)
				sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETWO_KIND;
			else if(line[1] == sys.exclude.EXC_ONETHREE)
				sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETHREE_KIND;
			else if(line[1] == sys.exclude.EXC_ONEFOUR)
				sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONEFOUR_KIND;
		}
		else if(line[0] == "Ewald")
		{
			sys.elect.ewald = checkBool(line[1]);
			sys.elect.readEwald = true;
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
		}
		else if(line[0] == "1-4scaling")
		{
			sys.elect.oneFourScale = stringtod(line[1]);
		}
		else if(line[0] == "Dielectric")
		{
			sys.elect.dielectric = stringtod(line[1]);
		}
		else if(line[0] == "RunSteps")
		{
			sys.step.total = stringtoi(line[1]);
		}
		else if(line[0] == "EqSteps")
		{
			sys.step.equil = stringtoi(line[1]);
		}
		else if(line[0] == "AdjSteps")
		{
			sys.step.adjustment = stringtoi(line[1]);
		}
		else if(line[0] == "DisFreq")
		{
			sys.moves.displace = stringtod(line[1]);
		}
		else if(line[0] == "IntraSwapFreq")
		{
			sys.moves.intraSwap = stringtod(line[1]);
		}
		else if(line[0] == "RotFreq")
		{
			sys.moves.rotate = stringtod(line[1]);
		} 
#ifdef VARIABLE_VOLUME
		else if(line[0] == "VolFreq")
		{
			sys.moves.volume = stringtod(line[1]);
		}
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
		else if(line[0] == "SwapFreq")
		{
#if ENSEMBLE == NVT
		        sys.moves.transfer = 0.000;
#else
			sys.moves.transfer = stringtod(line[1]);
#endif
		}
#endif
		else if(line[0] == "BoxDim") {
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
		else if(line[0] == "ChemPot") {
			if(line.size() != 3)
			{
				std::cout << "Error: Chemical potential parameters have not been specified!" << std::endl;
				exit(0);
			}
			std::string resName = line[1];
			double val = stringtod(line[2]);
			sys.chemPot.cp[resName] = val;
		} 
#endif
		else if(line[0] == "OutputName") 
		{
			out.statistics.settings.uniqueStr.val = line[1];
		} 
		else if(line[0] == "CoordinatesFreq") 
		{
			out.state.settings.enable = checkBool(line[1]);
			out.state.settings.frequency = stringtoi(line[2]);
		} 
		else if(line[0] == "RestartFreq") 
		{
			out.restart.settings.enable = checkBool(line[1]);
			out.restart.settings.frequency = stringtoi(line[2]);
		} 
		else if(line[0] == "ConsoleFreq") 
		{
			out.console.enable = checkBool(line[1]);
			out.console.frequency = stringtoi(line[2]);
		} 
		else if(line[0] == "BlockAverageFreq") 
		{
			out.statistics.settings.block.enable = checkBool(line[1]);
			out.statistics.settings.block.frequency = stringtoi(line[2]);
		} 
		else if(line[0] == "FluctuationFreq") 
		{
			out.statistics.settings.fluct.enable = checkBool(line[1]);
			out.statistics.settings.fluct.frequency = stringtoi(line[2]);
		} 
#if ENSEMBLE == GCMC
		else if(line[0] == "HistogramFreq") 
		{
			out.statistics.settings.hist.enable = checkBool(line[1]);
			out.statistics.settings.hist.frequency = stringtoi(line[2]);
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
		else if(line[0] == "OutEnergy") {
			out.statistics.vars.energy.block = checkBool(line[1]);
			out.statistics.vars.energy.fluct = checkBool(line[2]);
			out.statistics.vars.energy.hist = checkBool(line[3]);
		} else if(line[0] == "OutPressure") {
			out.statistics.vars.pressure.block = checkBool(line[1]);
			out.statistics.vars.pressure.fluct = checkBool(line[2]);
			out.statistics.vars.pressure.hist = checkBool(line[3]);
		}
#ifdef VARIABLE_PARTICLE_NUMBER
		else if(line[0] == "OutMolNum") {
			out.statistics.vars.molNum.block = checkBool(line[1]);
			out.statistics.vars.molNum.fluct = checkBool(line[2]);
			out.statistics.vars.molNum.hist = checkBool(line[3]);
		} else if(line[0] == "OutAcceptAngles") {
			out.statistics.vars.acceptAngles.block = checkBool(line[1]);
			out.statistics.vars.acceptAngles.fluct = checkBool(line[2]);
			out.statistics.vars.acceptAngles.hist = checkBool(line[3]);
		}
#endif
#ifdef VARIABLE_DENSITY
		else if(line[0] == "OutDensity") {
			out.statistics.vars.density.block = checkBool(line[1]);
			out.statistics.vars.density.fluct = checkBool(line[2]);
			out.statistics.vars.density.hist = checkBool(line[3]);
		}
#endif
#ifdef VARIABLE_VOLUME
		else if(line[0] == "OutVolume") {
			out.statistics.vars.volume.block = checkBool(line[1]);
			out.statistics.vars.volume.fluct = checkBool(line[2]);
			out.statistics.vars.volume.hist = checkBool(line[3]);
		}
#endif
		else{
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
	out.state.files.psf.name = prefix + out.statistics.settings.uniqueStr.val + ".psf"; /*ADDED PREFIX*/ 
	for(int i = 0; i<BOX_TOTAL; i++)
	{
		if(i==0)
			out.state.files.pdb.name[0] = prefix + out.statistics.settings.uniqueStr.val + "_BOX_0.pdb"; /*ADDED PREFIX*/ 
		else if(i==1)
			out.state.files.pdb.name[1] = prefix + out.statistics.settings.uniqueStr.val + "_BOX_1.pdb"; /*ADDED PREFIX*/ 
	}
	out.state.files.seed.name = prefix + out.statistics.settings.uniqueStr.val + ".dat"; /*ADDED PREFIX*/ 
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
	if(in.restart.enable == true && in.restart.step == ULONG_MAX)
	{
		std::cout << "Error: Restart step is needed!" << std::endl;
		exit(0);
	}
	if(in.restart.enable == false && in.restart.step != ULONG_MAX)
	{
		std::cout << "Warning: Restart step will not be used!" << std::endl;
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
	if(sys.elect.ewald && (sys.elect.tolerance == DBL_MAX))
	{
		std::cout << "Error: Tolerance has not been specified for Ewald summation!" << std::endl;
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
		std::cout << "Error: Temperature has not been defined!" << std::endl;
		exit(0);
	}
	if(out.statistics.settings.uniqueStr.val == "")
	{
		std::cout<< "Error: Output name has not been specified!" << std::endl;
		exit(0);
	}
	if(out.state.settings.enable && out.state.settings.frequency == ULONG_MAX)
	{
		std::cout << "Error: Coordinates output frequency has not been specified!" << std::endl;
		exit(0);
	}
	if(out.restart.settings.enable && out.restart.settings.frequency == ULONG_MAX)
	{
		std::cout << "Error: Restart coordinate output frequency has not been specified!" << std::endl;
		exit(0);
	}
	if(out.console.enable && out.console.frequency == ULONG_MAX)
	{
		std::cout << "Error: Console output frequency has not been specified!" << std::endl;
		exit(0);
	}
	if(out.statistics.settings.block.enable && out.statistics.settings.block.frequency == ULONG_MAX)
	{
		std::cout << "Error: Block average output frequency has not been specified!" << std::endl;
		exit(0);
	}
	if(out.statistics.settings.fluct.enable && out.statistics.settings.fluct.frequency == ULONG_MAX)
	{
		std::cout << "Error: Fluctuation output frequency has not been specified!" << std::endl;
		exit(0);
	}
	if(out.console.frequency == ULONG_MAX)
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
	if(out.restart.settings.frequency == ULONG_MAX)
	{
		out.restart.settings.frequency = (ulong)sys.step.total;
		std::cout << "Warning: By default restart coordinate output frequency has been set to " << out.restart.settings.frequency << "!" << std::endl;
	}
	if(out.state.settings.frequency == ULONG_MAX)
	{
		out.state.settings.frequency = (ulong)sys.step.total / 10;
		std::cout << "Warning: By default coordinate output frequency has been set to " << out.state.settings.frequency << "!" << std::endl;
	}
	if(out.statistics.settings.block.frequency == ULONG_MAX)
	{
	  
		out.statistics.settings.block.frequency = (ulong)sys.step.total / 100;
		std::cout << "Warning: By default block average output frequency has been set to " << out.statistics.settings.block.frequency << "!" << std::endl;
	}
	if(out.statistics.settings.fluct.frequency == ULONG_MAX)
	{
	  if(sys.step.total > 1000)
	  {
		out.statistics.settings.fluct.frequency = (ulong)sys.step.total / 1000;
	  }
	  else
	  {
	        out.statistics.settings.fluct.frequency = (ulong)sys.step.total / 100;
	  }
		std::cout << "Warning: By default fluctuation output frequency has been set to " << out.statistics.settings.fluct.frequency << "!" << std::endl;
	}
#if ENSEMBLE == GCMC
	if(out.statistics.settings.hist.enable && out.statistics.settings.hist.frequency == ULONG_MAX)
	{
		std::cout << "Error: Histogram output frequency has not been specified!" << std::endl;
		exit(0);
	}
	if(out.statistics.settings.hist.frequency == ULONG_MAX)
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
#ifdef VARIABLE_PARTICLE_NUMBER
	if(!out.statistics.settings.block.enable && out.statistics.vars.molNum.block)
	{
		std::cout<< "Warning: Output block average has been set off, molecule number output for block average will be ignored!" << std::endl;
	}
	if(!out.statistics.settings.block.enable && out.statistics.vars.acceptAngles.block)
	{
		std::cout<< "Warning: Output block average has been set off, accept angles output for block average will be ignored!" << std::endl;
	}
#endif
#ifdef VARIABLE_DENSITY
	if(!out.statistics.settings.block.enable && out.statistics.vars.density.block)
	{
		std::cout<< "Warning: Output block average has been set off, density output for block average will be ignored!" << std::endl;
	}
#endif
#ifdef VARIABLE_VOLUME
	if(!out.statistics.settings.block.enable && out.statistics.vars.volume.fluct)
	{
		std::cout<< "Warning: Output block average has been set off, volume ouput for block average will be ignored!" << std::endl;
	}
#endif
	if(!out.statistics.settings.fluct.enable && out.statistics.vars.energy.fluct)
	{
		std::cout<< "Warning: Output fluctuation has been set off, energy ouput for fluctuation will be ignored!" << std::endl;
	}
	if(!out.statistics.settings.fluct.enable && out.statistics.vars.pressure.fluct)
	{
		std::cout<< "Warning: Output fluctuation has been set off, pressure ouput for fluctuation will be ignored!" << std::endl;
	}
#ifdef VARIABLE_PARTICLE_NUMBER
	if(!out.statistics.settings.fluct.enable && out.statistics.vars.molNum.fluct)
	{
		std::cout<< "Warning: Output fluctuation has been set off, molecule number ouput for fluctuation will be ignored!" << std::endl;
	}
	if(!out.statistics.settings.fluct.enable && out.statistics.vars.acceptAngles.fluct)
	{
		std::cout<< "Warning: Output fluctuation has been set off, accept angles ouput for fluctuation will be ignored!" << std::endl;
	}
#endif
#ifdef VARIABLE_DENSITY
	if(!out.statistics.settings.fluct.enable && out.statistics.vars.density.fluct)
	{
		std::cout<< "Warning: Output fluctuation has been set off, density ouput for fluctuation will be ignored!" << std::endl;
	}
#endif
#ifdef VARIABLE_VOLUME
	if(!out.statistics.settings.fluct.enable && out.statistics.vars.volume.fluct)
	{
		std::cout<< "Warning: Output fluctuation has been set off, volume ouput for fluctuation will be ignored!" << std::endl;
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
