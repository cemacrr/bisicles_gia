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
// Transformation.H
// Abstract class defining a co-ordinate transformation, plus some concrete classes 
//===========================================================================
#ifndef _TRANSFORMATION_H_
#define _TRANSFORMATION_H_

#include "RealVect.H"
#include "NamespaceHeader.H"

/// abstract class for co-ordinate Transformations.
class Transformation
{
public:
 
  /// given x, compute x'
  virtual RealVect transform(const RealVect& x) const=0;
  
  virtual ~Transformation(){};
};



class PolarStereographicCartesianToPolarTransformation : public Transformation
{
  
/*
from: http://www.remotesensing.org/geotiff/proj_list/polar_stereographic.html

For the forward transformation from latitude and longitude:

E = FE + rho sin(lon - lon0)
N = FN - rho cos(lon - lon0)
where
rho = 2 a ko t /{[((1+e)^(1+e)) ((1-e)^(1-e))]^0.5}
t = tan (pi/4 - lat/2) / [(1-esin(lat) ) / (1 + e sin(lat))]^(e/2)

For the reverse transformation:

lat = chi+ (e^2/2 + 5e^4/24 + e^6/12 + 13e^8/360) sin(2 chi) 
+ (7e^4/48 + 29e^6/240 + 811e^8/11520) sin(4 chi)
+ (7e^6/120 +  81e^8/1120) sin(6 chi)  + (4279e^8/161280) sin(8 chi)

lon = lon0+ arctan [(E-FE) / (FN-N)]

where chi  = pi/2 - 2 arctan t
t   =  rho [((1+e)^(1+e)) ((1-e)^(1-e))]^0.5} / 2 a ko
rho = [(E-FE)^2  + (N - FN)^2]^0.5

ko == scale factor at natural origin (usually 1.0)
e = ellipsoid eccentricity (0.08181922 for earth)
a = equatorial radius (6.3781370e+6 m for earth)

*/ 


  Real m_e, m_a, m_k0, m_lon0;
  RealVect m_fx;

public:
  PolarStereographicCartesianToPolarTransformation
  (const Real& a_e, const Real& a_a, const Real& a_k0, const Real& a_lon0, 
   const RealVect& a_fx)
    :m_e(a_e),m_a(a_a),m_k0(a_k0),m_lon0(a_lon0),m_fx(a_fx)
  {}
  
  virtual RealVect transform(const RealVect& x) const
  {
    const Real& e = m_e;
    Real pi = M_PI;
    Real rho = std::sqrt( std::pow(x[0]-m_fx[0],2) +  std::pow(x[1]-m_fx[1],2));
    Real t = 0.5 * rho * std::sqrt ( std::pow(1.0+e,1.0+e) * std::pow(1.0-e,1.0-e)) /( 2.0 * m_a * m_k0);
    Real chi = 0.5*pi - 2.0*std::atan(t);
    Real lat = chi  
      + std::sin(2.0*chi) * (0.5 * std::pow(e,2) + 5.0*std::pow(e,4)/24.0 +  std::pow(e,6)/12.0 + 13.0*std::pow(e,8)/360.0)
      + std::sin(4.0*chi) * (7.0*std::pow(e,4)/48.0 +  29.0*std::pow(e,6)/240.0 + 811.0*std::pow(e,8)/11520.0)
      + std::sin(6.0*chi) * (7.0*std::pow(e,6)/120.0 + 81.0*std::pow(e,8)/1120.0)
      + std::sin(8.0*chi) * (4729.0*std::pow(e,8)/161280.0);

   
    Real lon =  m_lon0 + std::atan2 ( x[0]-m_fx[0] , m_fx[1]-x[1] );
    
    RealVect r ( x) ;
    r[0] = lat;
    r[1] = lon;

  return r;

  }

  
};


#include "NamespaceFooter.H"
#endif
