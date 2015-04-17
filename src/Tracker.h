/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef TRACKER_H
#define TRACKER_H

#include <limits>
#include <map>
#include <cmath>

#include "Globals.h"
#include "InputFile.h"

#if ENSEMBLE == GEMC
#define NUM_TRKD_BOXES 2
#else
#define NUM_TRKD_BOXES 1
#endif

template <typename T>
struct SvEntry
{
   unsigned long step;
   T val;
};

#define BIN_ARRAY true
#define BIN_MAP false

template <typename T>
struct Tracker
{  
 public:

   Tracker(InputFile & conf, std::string const& vName, 
	   std::string const& ivStr, const T defVal,
	   const unsigned long totSteps, const T rangeMax,
	   T const*const rangeMin =NULL) : 
   currStep(0), varName(vName), indValStr(ivStr), defaultVal(defVal),
      binKind(BIN_ARRAY), min(NULL), max(NULL)
   {     
      //Read blk vars
      BlkInit(conf,totSteps);

      //Read fluct vars;
      FluctInit(conf);

      //Read Histogram variables
      HistInit(conf, rangeMin, &rangeMax);
   }

   Tracker(InputFile & conf, std::string const& vName, 
	   std::string const& ivStr, const T defVal,
	   const unsigned long totSteps) : 
   currStep(0), varName(vName), indValStr(ivStr), defaultVal(defVal), 
      binKind(BIN_MAP), min(NULL), max(NULL)
   {     
      //Read blk vars
      BlkInit(conf,totSteps);

      //Read fluct vars;
      FluctInit(conf);

      //Read Histogram variables
      HistInit(conf, NULL, NULL);
   }

   ~Tracker()
   {
      //Do final check for writes
      if (blkOn)
	 BlkCheckForWrt(true);
      if (fluctOn)
	 FluctCheckForWrt(true);
      if (histOn)
	 HistCheckForWrt(true);

      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
      {
	 if (blkOn)
	    delete[] blkAvg[b];
	 if (binKind==BIN_ARRAY&&histOn)
	    delete[] histArr[b];
	 if (fluctOn)
	    delete[] entries[b];
      }
      if (histOn)
      {
	 delete[] min;
	 delete[] max;
      }
   }

   inline bool SampleOK(const T chk)
   { return true; }

   inline void Sample(T const*const val, const unsigned long step,
		      const bool suppressWrt = false)
   {
      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	 if (!SampleOK(val[b])) 
	    return;
      if (step != currStep)
	 currStep = step;
      if (blkOn)
	 BlkSample(val,suppressWrt);
      if (fluctOn)
	 FluctSample(val);
      if (histOn)
	 HistSample(val,suppressWrt);
   }

   //Sample from only one box.
   inline void Sample(const T val, const unsigned int b,
		      const unsigned long step,
		      const bool suppressWrt = false)
   {
      //In GCMC don't track second box.
      if (b >= NUM_TRKD_BOXES || !SampleOK(val) )
	 return;
      if (step != currStep)
	 currStep = step;
      if (blkOn)
	 BlkSample(val,b,suppressWrt);
      if (fluctOn)
	 FluctSample(val,b);
      if (histOn)
	 HistSample(val,b,suppressWrt);
   }

   inline void CheckForWrt(const unsigned long step)
   {
      if (step != currStep)
	 currStep = step;
      if (blkOn)
	 BlkCheckForWrt();
      if (fluctOn)
	 FluctCheckForWrt();
      if (histOn)
	 HistCheckForWrt();
   }

 private:
   ///////////////////
   //OVERALL
   ////////////////////

   //Number of samples collected
   unsigned long currStep;
   //Name of tracked var;
   std::string varName;
   //String containing independent values for current sim.
   std::string indValStr;
   //Default value to set new bins/etc. to
   T defaultVal;

   inline std::string GetFileName(std::string const& knd)
   {
      std::string fNm;
      fNm = knd;
      fNm += "_";
      fNm += varName;
      fNm += "_";
      fNm += indValStr;
      fNm += ".dat";
      return fNm;
   }

   inline void GetFileNames(std::string * fNm, std::string const& knd)
   {
      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
      {
	 fNm[b] = knd;
	 fNm[b] += "_";
	 fNm[b] += varName;
	 fNm[b] += "_";
	 fNm[b] += indValStr;
	 fNm[b] += "_Box_";
	 fNm[b] += str<unsigned int>(b);
	 fNm[b] += ".dat";
      }
   }

   ///////////////////////
   //BLOCK Avg. Vars.
   ///////////////////////

   //Whether to use block averages
   bool blkOn;
   //How fine to divide into blocks.
   unsigned int stepsPerBlk;
   //Block average totals
   T * blkAvg[NUM_TRKD_BOXES];
   //Running total, will be written to current blk either at end of blk
   //or however often is needed to prevent overflow.
   T currBlkTemp[NUM_TRKD_BOXES];
   //Running total of # of samples to divide by.
   unsigned samplesSinceAdd[NUM_TRKD_BOXES];
   //Name of file to write to
   std::string blkFileNm;
   //If this has a value...
   bool blkHasVal[NUM_TRKD_BOXES];
   //Whether to overwrite blk files.
   bool blkFirstWrt;

   inline void DoSubAvg(const unsigned int b)
   {
      unsigned int idx = currStep/stepsPerBlk;   
      blkAvg[b][idx] += currBlkTemp[b]/samplesSinceAdd[b];
      currBlkTemp[b] = 0;
      samplesSinceAdd[b] = 0;
   }

   inline void BlkCheckForWrt(bool forceWrt=false)
   {
      //If last sample of block, write to file.
      if (((currStep+1)%stepsPerBlk==0) != forceWrt)
      {
	 unsigned int idx = currStep/stepsPerBlk;
	 //Write out
	 std::ofstream writer;
	 writer.open(blkFileNm.c_str(),
		     ((blkFirstWrt)?std::ios::out:std::ios::app));
      
	 writer << currStep;
	 for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	 {
	    if (samplesSinceAdd[b]>0)
	       DoSubAvg(b);
	    WriteSpace(writer);
	    if(blkHasVal[b])
	       Write<T>(writer,blkAvg[b][idx],12, std::ios::left, 12); 
	    else
	       Write<T>(writer,defaultVal,12, std::ios::left, 12); 
	 }
	 WriteEndl(writer);
	 writer.close();
	 for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	    blkHasVal[b] = false;
	 blkFirstWrt=false;
      }
   }

   inline void BlkInit(InputFile & conf, const unsigned long totSteps)
   {
      conf.Goto("BlkDump",varName);
      blkOn = conf.GetBool();
      if (blkOn)
      {
	 blkFirstWrt = true;
	 blkFileNm = GetFileName("Blk");
	 stepsPerBlk = conf.GetValue<unsigned int>();
	 for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	 {
	    blkHasVal[b] = false;
	    currBlkTemp[b] = 0;
	    samplesSinceAdd[b] = 0;
	    blkAvg[b] = new T[totSteps/stepsPerBlk+1];
	    for (unsigned int bl = 0; bl < totSteps/stepsPerBlk; bl++)
	       blkAvg[b][bl] = defaultVal;
	 }
      }
   }
   inline void BlkSample(const T val, const unsigned int b,
			 const bool suppressWrt = false)
   {
      currBlkTemp[b] += val;
      samplesSinceAdd[b]++;
      blkHasVal[b] = true;
	 
      //If overflow possible, store results and zero out intermediates
      if (std::numeric_limits<T>::max()/2 < currBlkTemp[b])
	 DoSubAvg(b);

      if (!suppressWrt)
	 BlkCheckForWrt();
   }

   inline void BlkSample(T const*const val, const bool suppressWrt = false)
   {
      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
      {
	 currBlkTemp[b] += val[b];
	 samplesSinceAdd[b]++;
	 blkHasVal[b] = true;
	 //If overflow possible, store results and zero out intermediates
	 if (std::numeric_limits<T>::max()/2 < currBlkTemp[b])
	    DoSubAvg(b);
      }
      if (!suppressWrt)
	 BlkCheckForWrt();
   }

   /////////////////////////
   //FLUCTUATION Vars.
   /////////////////////////

   //Whether to use fluctuation file
   bool fluctOn;
   //Set of saved samples.
   SvEntry<T> * entries[NUM_TRKD_BOXES];
   //Number of samples
   unsigned int samplesPerEntry;
   //How many samples to store before write.
   unsigned int entriesPerWrt;
   //Name of file to write to
   std::string fluctFileNm[NUM_TRKD_BOXES];
   //If this has a value...
   unsigned int sampleCnt[NUM_TRKD_BOXES];
   //If this has a value...
   unsigned int entriesCnt[NUM_TRKD_BOXES];
   //Whether to overwrite blk files.
   bool fluctFirstWrt[NUM_TRKD_BOXES];

   inline void FluctCheckForWrt(bool forceWrt=false)
   {
      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
      {
	 if ((entriesCnt[b] == entriesPerWrt) != forceWrt)
	 {
	    if (varName.compare("Volume")==0)
	       printf("hi!\n");
	    std::ofstream writer;
	    writer.open(fluctFileNm[b].c_str(),
			(fluctFirstWrt[b]?std::ios::out:std::ios::app));
	    for (unsigned int e = 0; e<entriesCnt[b] ; e++)
	       writer << entries[b][e].step << " " << entries[b][e].val 
		      << std::endl;
	    writer.close();
	    sampleCnt[b]=0;
	    entriesCnt[b]=0;
	    fluctFirstWrt[b] = false;
	 }
      }
   }

   inline void FluctInit(InputFile & conf)
   {
      conf.Goto("FluctDump",varName);
      fluctOn = conf.GetBool();
      if (fluctOn)
      {
         samplesPerEntry = conf.GetValue<unsigned int>();
         entriesPerWrt = conf.GetValue<unsigned int>();
	 for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	 {
	    entries[b] = new SvEntry<T>[entriesPerWrt];
	    GetFileNames(fluctFileNm,"Fluct");
	    sampleCnt[b] = 0;
	    entriesCnt[b] = 0;
	    fluctFirstWrt[b] = true;
	 }
      }
   }

   inline void FluctSample(const T val, const unsigned int b)
   {
      sampleCnt[b]++;
      
      if (sampleCnt[b]%samplesPerEntry==0)
      {
	 entries[b][entriesCnt[b]].step = currStep;
	 entries[b][entriesCnt[b]].val = val;
	 entriesCnt[b]++;
	 sampleCnt[b] = 0;
      }
      
      FluctCheckForWrt();
   }

   inline void FluctSample(T const*const val)
   {
      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	 FluctSample(val[b],b);
   }


   /////////////////////
   //HISTOGRAM Vars.
   ////////////////////
   
   //Whether to use histogram
   bool histOn;
   //Flag to control which bin is used;
   bool binKind;
   //Type 1: 
   //Known range, known # bins
   //... or ...
   //Known range, known bin size
   unsigned int * histArr [NUM_TRKD_BOXES];
   unsigned int binCnt [NUM_TRKD_BOXES];
   //Type 2:
   //Known bin size only
   std::map<T,unsigned int> histMap[NUM_TRKD_BOXES];
   T * min, * max;
   //Size of our bin.
   T binSz;
   //Steps per histogram dump
   unsigned int stepsPerHist;
   //Name of file to write to
   std::string histFileNm[NUM_TRKD_BOXES];   
   //If this has a value...
   bool histHasVal[NUM_TRKD_BOXES];

   inline void HistCheckForWrt(bool forceWrt=false)
   {
      T val = defaultVal;
      std::ofstream writer;
      if ( ((currStep+1)%stepsPerBlk==0&&binCnt!=NULL) != forceWrt )
      {
	 for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	 {
	    if (histHasVal[b])
	    {
	       writer.open(histFileNm[b].c_str(),std::ios::out);
	       if (binKind==BIN_ARRAY)
	       {
		  val = min[b];
		  for (unsigned int bin = 0; bin < binCnt[b]; bin++)
		  {
		     writer << val << " " << histArr[b][bin] << std::endl;
		     val += binSz;
		  }
	       }
	       else
	       {
		  for(typename std::map< T, unsigned int >::iterator it = 
			 histMap[b].begin(); 
		      it != histMap[b].end();
		      it++)
		     writer << it->first << " " << it->second << std::endl;
	       }
	       writer.close();
	    }
	 }
      }
   }

   inline void HistInit(InputFile & conf, T const*const rangeMin,
			T const*const rangeMax)
   {
      conf.Goto("HistDump",varName);
      histOn = conf.GetBool();
      if (histOn)
      {
	 min = new T[NUM_TRKD_BOXES];
	 max = new T[NUM_TRKD_BOXES];
	 for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	       histHasVal[b] = false;
	 GetFileNames(histFileNm,"Hist");
	 stepsPerHist = conf.GetValue<unsigned int>();
	 binSz = conf.GetValue<T>();
	 if (binKind==BIN_ARRAY)
	 {
	    for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	    {
	       min[b] = (rangeMin!=NULL?*rangeMin:defaultVal); 
	       max[b] = *rangeMax;
	       binCnt[b] = max[b]/binSz;
	       histArr[b] = new unsigned int[binCnt[b]];
	       for (unsigned int bin = 0; bin < binCnt[b]; bin++)
		  histArr[b][bin] = 0;
	    }
	 }
	 else
	    for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
	       binCnt[b] = 0;
      }
   }

   inline void HistSample(const T val, const unsigned int b,
			  const bool suppressWrt)
   {
      histHasVal[b] = true;
      if (binKind==BIN_MAP)
      {
	 //Truncate to binSz.
	 T binTgt = ((T)((int)(val/binSz)))*binSz;
	 if (binCnt[b]==0)
	 {
	    if (min==NULL)
	    {
	       min = new T[NUM_TRKD_BOXES];
	       max = new T[NUM_TRKD_BOXES];
	    }
	    min[b] = max[b] = binTgt;
	    histMap[b].insert(std::pair<T,unsigned int>(binTgt,1));
	    binCnt[b]++;
	 }
	 else if (histMap[b].find(binTgt)==histMap[b].end())
	 {
	    histMap[b].insert(std::pair<T,unsigned int>(binTgt,1));
	    binCnt[b]++;
	    min[b] = (binTgt < min[b]?binTgt:min[b]);
	    max[b] = (binTgt > max[b]?binTgt:max[b]);
	 }
	 else
	 {
	    histMap[b][binTgt]++;
	 }
      }
      else
      {	 
	 unsigned int idx = (unsigned int)((val-min[b])/binSz);
	 histArr[b][idx]++;
      }
      
      if (!suppressWrt)
	 HistCheckForWrt();
   }

   inline void HistSample(T const*const val, const bool suppressWrt)
   {
      for (unsigned int b = 0; b < NUM_TRKD_BOXES; b++)
      {
	 histHasVal[b] = true;
	 if (binKind==BIN_MAP)
	 {
	    //Truncate to binSz.
	    T binTgt = ((T)((int)(val[b]/binSz)))*binSz;
	    if (binCnt[b]==0)
	    {
	       if (min==NULL)
	       {
		  min = new T[NUM_TRKD_BOXES];
		  max = new T[NUM_TRKD_BOXES];
	       }
	       min[b] = max[b] = binTgt;
	       histMap[b].insert(std::pair<T,unsigned int>(binTgt,1));
	       binCnt[b]++;
	    }
	    else if (histMap[b].find(binTgt)==histMap[b].end())
	    {
	       histMap[b].insert(std::pair<T,unsigned int>(binTgt,1));
	       binCnt[b]++;
	       min[b] = (binTgt < min[b]?binTgt:min[b]);
	       max[b] = (binTgt > max[b]?binTgt:max[b]);
	    }
	    else
	    {
	       histMap[b][binTgt]++;
	    }
	 }
	 else
	 {	 
	    unsigned int idx = (unsigned int)((val[b]-min[b])/binSz);
	    histArr[b][idx]++;
	 }
      }
      if (!suppressWrt)
	 HistCheckForWrt();
   }
};
   
template <>
inline bool Tracker<double>::SampleOK(const double chk)
{
   return std::isfinite(chk)!=0 &&
      chk > (std::numeric_limits<double>::max()/10.0*-1.0) &&
      chk < (std::numeric_limits<double>::max()/10.0);
}

#endif /*TRACKER_H*/

