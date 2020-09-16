#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MONOMIALIF_H_
#define _MONOMIALIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"
#include "OneDIntegrator.H"

#include "NamespaceHeader.H"

///
/**
   class to implement 1-D function that gets integrated 
   to compute EBMoments for monomials
*/
class MonomialEBMomentFunc : public BaseOneDFunc
{
public:
  ///
  /**
       Q = power of monomial
       p = moment being calculated
   */
  MonomialEBMomentFunc(const IvSDMinOne& a_Q,
                       const IvSpaceDim& a_p)
  {
    m_Q = a_Q;
    m_p = a_p;
  }


  virtual ~MonomialEBMomentFunc()
  {; }

  virtual Real value(const Real& location) const;

private:
  //power of monomial
  IvSDMinOne m_Q;
  //moment being calculated
  IvSpaceDim m_p;
  //strong construction is best when it can be used
  MonomialEBMomentFunc()
  { ;  }
  
};
///
/**
   This implicit function specifies a Monomial.
*/
class MonomialIF: public BaseIF
{
public:
  ///
  /**
     in 2d
     y = (x)^p[0] 
     in 3d
     z = [(x)^p[0]*(y)^p[1]]

     p == a_power
  */
  MonomialIF(const IndexTM<int, SpaceDim-1> & a_power,
             const bool                     & a_inside);

  /// Destructor
  virtual ~MonomialIF()
  {
  }

  virtual Real value(const IndexTM<int,SpaceDim> & a_partialDerivative,
                     const IndexTM<Real,SpaceDim>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,SpaceDim>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  ///return int x^p dV for the vof
  virtual IndMomSpaceDim getExactVolumeMoments(const VolIndex    & a_vof,
                                               const Real        & a_dx) const;
                                          
  
  ///return int_eb x^p dA for the eb
  virtual IndMomSpaceDim getExactEBMoments(const VolIndex   & a_vof,
                                           const Real       & a_dx) const;
                                          
                                          
  ///return int_eb x^p n_i dA for the eb
  virtual IndMomSpaceDim getExactEBNormalMoments(const VolIndex   & a_vof,
                                                 const Real       & a_dx,
                                                 const int        & a_ni) const;
                                          
  ///
  virtual IndMomSpaceDim getExactEBNormalPartialDerivs(const VolIndex   & a_vof,
                                                       const Real       & a_dx,
                                                       const int        & a_ni) const;

  ///
  Real
  getZerothDerivativeOfXNormal(const Real& a_xloc) const;
  ///
  Real
  getFirstDerivativeOfXNormal(const Real& a_xloc) const;
  ///
  Real
  getSecondDerivativeOfXNormal(const Real& a_xloc) const;
  ///
  Real
  getThirdDerivativeOfXNormal(const Real& a_xloc) const;
  ///
  Real
  getFourthDerivativeOfXNormal(const Real& a_xloc) const;
  ///
  Real
  getZerothDerivativeOfYNormal(const Real& a_xloc) const; 
  ///
  Real
  getFirstDerivativeOfYNormal(const Real& a_xloc) const;
  ///
  Real
  getSecondDerivativeOfYNormal(const Real& a_xloc) const;
  ///
  Real
  getThirdDerivativeOfYNormal(const Real& a_xloc) const;
  ///
  Real
  getFourthDerivativeOfYNormal(const Real& a_xloc) const;


  ///return int x^p  dA for the face
  virtual IndMomSDMinOne getExactFaceMoments(const FaceIndex  & a_face,
                                             const Real       & a_dx) const;
                                          
  void getTloThi(bool& a_noEB,
                 Real& a_tlo, 
                 Real& a_thi,
                 const VolIndex   & a_vof,
                 const Real       & a_dx) const;

  Real getEBNormMomX(const Real        & a_tlo,
                     const Real        & a_thi,
                     const IvSpaceDim  & a_p) const;

  Real getEBNormMomY(const Real        & a_tlo,
                     const Real        & a_thi,
                     const IvSpaceDim  & a_p) const;

  Real xToTheQ(const IndexTM<Real,SpaceDim>& a_point) const;

  IndMomSpaceDim 
  getUncutVolumeSection(const RealVect& a_xlo, 
                        const RealVect& a_xhi, 
                        const RealVect& a_xbar,
                          const Real    & a_dx) const;

  IndMomSpaceDim 
  getCutVolumeSection(const RealVect & a_xlo, 
                      const RealVect & a_xhi, 
                      const RealVect & a_xbar,
                      const Real     & a_dx) const;

  void shiftToXBar(IndMomSpaceDim& a_moment,
                   const RealVect& a_xbar) const;
protected:

  IndexTM<int, SpaceDim-1>  m_power;
  bool                      m_inside;

private:
  MonomialIF();

};

#include "NamespaceFooter.H"
#endif
