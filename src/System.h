#ifndef SYSTEM_H
#define SYSTEM_H

#include "EnsemblePreprocessor.h" //For VARIABLE_<QUANTITY> conditional defines
#include "CalculateEnergy.h" 
#include "Ewald.h" 

//Member variables
#include "EnergyTypes.h"
#include "Coordinates.h"
#include "PRNG.h"
#include "BoxDimensions.h"
#include "MoleculeLookup.h"
#include "MoveSettings.h"
#ifdef CELL_LIST
#include "CellList.h"
#endif

//Initialization variables
class Setup;
class StaticVals;
class MoveBase;

class System
{
 public:
   explicit System(StaticVals& statics);

   void Init(Setup const& setupData);

   //Runs move, picked at random
   void ChooseAndRunMove(const uint step);

   // return ewald
   Ewald * GetEwald()
   {
     return calcEwald;
   }

   //NOTE:
   //This must also come first... as subsequent values depend on obj.
   //That may be in here, i.e. Box Dimensions
   StaticVals & statV;

   //NOTE:
   //Important! These must come first, as other objects may depend
   //on their val for init!
   //Only include these variables if they vary for this ensemble...
#ifdef VARIABLE_VOLUME
   BoxDimensions boxDimensions;
#endif
#ifdef  VARIABLE_PARTICLE_NUMBER
   MoleculeLookup molLookup;
#endif

   //Use as we don't know where they are...
   BoxDimensions & boxDimRef;
   MoleculeLookup & molLookupRef;

   MoveSettings moveSettings;
   SystemPotential potential;
   Coordinates coordinates;
   COM com;

   CalculateEnergy calcEnergy;
   Ewald  *calcEwald;
   //Ewald ewaldEnergy;
   //Ewald & calcEwald;
#ifdef CELL_LIST
   CellList cellList;
#endif
   PRNG prng;

   //Procedure to run once move is picked... can also be called directly for
   //debugging...
   void RunMove(uint majKind, double draw, const uint step);

   ~System();

 private:
   void InitMoves();
   void PickMove(uint & kind, double & draw);
   uint SetParams(const uint kind, const double draw);
   uint Transform(const uint kind);
   void CalcEn(const uint kind);
   void Accept(const uint kind, const uint rejectState, const uint step);

   MoveBase * moves[mv::MOVE_KINDS_TOTAL];
};

#endif /*SYSTEM_H*/
