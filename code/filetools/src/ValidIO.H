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
// ValidIO.H
// Operations related to reading and writing data from/to valid regions
// of AMR Hierarchies
//===========================================================================
#ifndef _VALIDIO_H_
#define _VALIDIO_H_

#include "RealVect.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "netcdf.h"
#include "Transformation.H"
#include "NamespaceHeader.H"

class ValidData
{

  int m_nComp;
  Vector<RealVect> m_dx;
  RealVect m_x0;
  int m_nLevel;
  Vector<int> m_ratio;
  Vector<int> m_level;
  Vector<int> m_boxID;
  Vector<Vector<int> > m_iv;
  Vector<Vector< Real > > m_field;
  Vector<Vector<Real > > m_x;

  bool m_isDefined;

public:
  ValidData(int a_nComp,  const RealVect& a_crseDx, const Vector<int>& a_ratio, const RealVect& a_x0 = RealVect::Zero );
  
  void define(int a_nComp,  const RealVect& a_crseDx, const Vector<int>& a_ratio, const RealVect& a_x0 = RealVect::Zero );

  ValidData():m_isDefined(false){}

  void append(int a_lev, int a_boxID, const IntVect& a_iv, 
	      const Vector<Real>& a_data);
  
  const int nComp() const
  {
    return m_nComp;
  }

  const size_t nCell() const
  {
    return m_level.size();
  }

  const size_t nLevel() const
  {
    return m_nLevel;
  }

  const Vector<int>& ratio() const
  {
    return m_ratio;
  }
 
  const Vector<RealVect>& dx() const
  {
    return m_dx;
  }


  const Vector<int>& level() const
  {
    return m_level;
  }

  const Vector<int>& boxID() const
  {
    return m_boxID;
  }

  const Vector<int>& iv(int a_dir) const
  {
    return m_iv[a_dir];
  }

  const Vector<Real>& field(int a_comp) const
  {
    return m_field[a_comp];
  }
  const Vector<Real>& x(int a_dir) const
  {
    return m_x[a_dir];
  }

  ///fill  a_nodeCoord with the a_dir co-ordinate of node a_node. 
  /**
     \param Vector<Real>& a_nodeCoord will be resized to this.size()
     \param a_node[dir] = 0 implies lo side, 1 implies high side
  */
  void computeNodeCoord(Vector<Real>& a_nodeCoord, int a_dir, const IntVect& a_node) const;
    
};

namespace ValidIO 
{

  /// populate a ValidData from a  Vector<LevelData<FArrayBox>*>
  void BStoValid ( ValidData & , 
		   const Vector<LevelData<FArrayBox>*>&, 
		 const Vector<int>& a_ratio );

  ///write ValidData to a NetCDF-CF file
  void writeCF ( const std::string&,  const ValidData& , 
		 const Vector<std::string>&, 
		 const Transformation&);

  ///read ValidData from a NetCDF-CF file
  void readCF ( ValidData&, Vector<std::string>&, 
		 const std::string&);



}

#include "NamespaceFooter.H"
#endif
