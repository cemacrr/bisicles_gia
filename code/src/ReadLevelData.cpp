#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// ReadLevelData.cpp
//===========================================================================


#include "ReadLevelData.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "NamespaceHeader.H"

//fill a LevelData<FArrayBox> from data stored in an AMR file
//for each name in a_names look for either 
//(a) a single component LevelData named "name"
//(b) a multi-component LevelData, stored consecutively and beginning with name
void readLevelData(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_data,
		   Real& a_dx,
		   const std::string a_file,
		   const Vector<std::string>& a_names, 
		   int a_nComp)
{

  Vector<LevelData<FArrayBox>* > vectData;
  Vector<DisjointBoxLayout> vectGrids;
  Vector<int> vectRatio;
  Vector<std::string> names;
  Real dt = 0.0,time = 0.0;
  Box domBox;
  int numLevels;

  pout() << " attempting to open file  " << a_file << std::endl;

  int status = ReadAMRHierarchyHDF5
    (a_file,vectGrids,vectData,names,domBox,a_dx,dt,time,
     vectRatio,numLevels);

 CH_assert(status == 0);
 if (status != 0)
   MayDay::Error("failed to read file");
 CH_assert(vectData.size() == 1);
 if (vectData.size() != 1)
   MayDay::Error("bad data");
 
 
 Vector<Box> boxes;
 int max_box_size = 64;
 int block_factor = 8;
 domainSplit(domBox, boxes, max_box_size, block_factor);
 Vector<int> procAssign(boxes.size());
 LoadBalance(procAssign,boxes);

 DisjointBoxLayout grids(boxes, procAssign, domBox);

 for (int i =0; i < a_data.size(); ++i)
   a_data[i]->define(grids,a_nComp,IntVect::Zero);

 int read = 0;
 
 for (int j = 0; j < a_names.size(); j++)
   {
     pout() << " looking for variable " << a_names[j] << std::endl;
     for (int i = 0; i < names.size(); i++)
       {
	 if (names[i] == a_names[j])
	   {
	     pout() << " found variable " << names[i] << std::endl;
	     vectData[0]->copyTo(Interval(i,i+a_nComp-1),*a_data[j],Interval(0,a_nComp-1));
	     read++;
	   }
       }
   }

 CH_assert(read == a_names.size());
 if (!read)
   MayDay::Error("no data in AMR Hierarchy");

 for (int i = 0; i < vectData.size(); i++)
   {
     if (vectData[i] != NULL)
       {
	 delete vectData[i]; vectData[i] = NULL;
       }
   }

}


#include "NamespaceFooter.H"
