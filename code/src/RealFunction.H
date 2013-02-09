#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


/// RealFunction.H
/// ==============
/// defines virtual base template for real functions
#ifndef _REALFUNCTION_H_
#define _REALFUNCTION_H_
#include "CH_HDF5.H"
#include "NamespaceHeader.H"

// virtual base template which defines a Real function
template <typename T> 
class RealFunction
{
public:
  virtual ~RealFunction()
  {;}
  virtual Real operator()(const T&)=0;
};

// simple constant
template <typename T>
class ConstantRealFunction : public RealFunction<T>
{
  Real m_value;
public:
  ConstantRealFunction(Real a_value) : m_value(a_value)
  {;}
  Real operator()(const T&){return m_value;}
};

//linear function of a realVect
class LinearRealFunction : public RealFunction<RealVect>
{
  Real m_originValue;
  RealVect m_gradient;

public:

  Real operator()(const RealVect& x)
  {
    return m_originValue + (m_gradient * x).sum();
  }

};


class InclinedPlaneFunction : public RealFunction<RealVect>
{
  Real m_origin;
  RealVect m_slope;
public:
  InclinedPlaneFunction(Real a_origin, const RealVect& a_slope)
    : m_origin(a_origin),m_slope(a_slope)
  {;}

  Real operator()(const RealVect& x)
  {
     Real r =  m_origin + 
       D_TERM(x[0]*m_slope[0],+x[1]*m_slope[1],+x[2]*m_slope[2]);

     return r;

  }

};


class SymmetricInclinedPlaneFunction : public RealFunction<RealVect>
{
  Real m_origin;
  RealVect m_slope;
  RealVect m_symmetryPoint;
public:
  SymmetricInclinedPlaneFunction(Real a_origin, const RealVect& a_slope, 
                        const RealVect& a_symmetryPoint)
    : m_origin(a_origin),m_slope(a_slope), m_symmetryPoint(a_symmetryPoint)
  {;}

  Real operator()(const RealVect& x)
  {
    RealVect local_x(x);
    local_x -= m_symmetryPoint;
    D_TERM(local_x[0] = abs(local_x[0]);,
           local_x[1] = abs(local_x[1]);,
           local_x[2] = abs(local_x[2]);)


     Real r =  m_origin + 
       D_TERM(local_x[0]*m_slope[0],+x[1]*m_slope[1],+x[2]*m_slope[2]);

     return r;

  }

};


class GaussianFunction : public RealFunction<RealVect>
{
  RealVect m_center;
  RealVect m_radius;
  Real m_mag;
public:
  GaussianFunction(RealVect a_center, const RealVect& a_radius, Real a_mag)
    : m_center(a_center),m_radius(a_radius), m_mag(a_mag)
  {;}

  Real operator()(const RealVect& x)
  {
    RealVect loc(x);
    loc -= m_center;
    loc /= m_radius;
    Real radSqr = D_TERM(loc[0]*loc[0],+loc[1]*loc[1],+loc[2]*loc[2]);
    
    Real r =  m_mag*exp(-radSqr);
    
    return r;

  }

};

/// 1D Gaussian hump
class GaussianFunctionX : public RealFunction<RealVect>
{
  Real m_center;
  Real m_radius;
  Real m_mag;
public:
  GaussianFunctionX(Real a_center, const Real& a_radius, Real a_mag)
    : m_center(a_center),m_radius(a_radius), m_mag(a_mag)
  {;}

  Real operator()(const RealVect& x)
  {
    RealVect loc(x);
    loc -= m_center;
    loc /= m_radius;
    Real radSqr = loc[0]*loc[0];
    
    Real r =  m_mag*exp(-radSqr);
    
    return r;

  }

};


/// allow flowline data to be read from an hdf5 file
/// Flowline is parallel to the x axis
/// Data must contain n values f_i, where f_i = f(dx * (i + 1/2))
/// and the problem domain has 0 < x < (n - 1/2)*dx

class ExtrudedPieceWiseLinearFlowline : public RealFunction<RealVect>
{

  Vector<Real> m_flowlineF; 
  Real m_dx;
public:

  ExtrudedPieceWiseLinearFlowline(std::string a_file,  std::string a_set, Real a_dx)
    :m_dx(a_dx)
  {
#ifdef CH_USE_HDF5
    herr_t status;
    hid_t file_access = 0;
#ifdef CH_MPI
    file_access = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(file_access,  Chombo_MPI::comm, MPI_INFO_NULL);
#else
    file_access = H5P_DEFAULT;
#endif
    hid_t file_id = H5Fopen(a_file.c_str(), H5F_ACC_RDONLY, file_access);
    hid_t dset_id = H5Dopen(file_id, a_set.c_str());
    hid_t space_id = H5Dget_space(dset_id);
    CH_assert( H5Sis_simple(space_id) > 0);
    CH_assert( H5Sget_simple_extent_ndims(space_id) == 1);
    hsize_t n,nmax;
    status = H5Sget_simple_extent_dims(space_id,&n,&nmax);
    m_flowlineF.resize(n);
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&m_flowlineF[0]);
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
    status = H5Dclose(dset_id);
    status = H5Sclose(space_id);
    status = H5Fclose(file_id); 
#else
    MayDay::Error("ExtrudedPieceWiseLinearFlowline(std::string a_file,  std::string a_set, Real a_dx) requires CH_USE_HDF5");
#endif    
      } 

  Real operator()(const RealVect& a_x)
  {
    Real x = abs(a_x[0]);
    int i = int( x / m_dx);
  
    if (i + 1 < m_flowlineF.size())
      {
	Real w = (x-(0.5+Real(i))*m_dx)/m_dx;
	return m_flowlineF[i]*(1.0-w) + m_flowlineF[i+1]*w;
      }
    else
      {
	return m_flowlineF[m_flowlineF.size()-1];
      }
  }
  

};




#include "NamespaceFooter.H"
#endif
 
