#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MonomialIF.H"
#include "Factorial.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

///
Real 
sumAllDivisions(const Real* a_numers,
                const Real* a_denoms,
                const int   a_numSum) 
{
  Real retval = 0;
  for(int isum = 0; isum < a_numSum; isum++)
    {
      retval += a_numers[isum]/a_denoms[isum];
    }
  return retval;
}
/*******/
Real 
MonomialEBMomentFunc::
value(const Real& a_loc) const
{
  CH_assert(SpaceDim == 2);
  CH_assert(m_Q[0] >= 1);
  Real Px = m_p[0];
  Real Py = m_p[1];
  Real Q = m_Q[0];
  Real cpow = Max(2*Q-2., 0.);
  Real absCPrim = sqrt(1.0 + Q*Q*pow(a_loc, cpow));
  Real retval   = absCPrim*pow(a_loc, Px + Q*Py);
  return retval;
}

/*******/
MonomialIF::
MonomialIF(const IndexTM<int, SpaceDim-1> &     a_power,
           const bool                     &     a_inside)
  :BaseIF()
{
  //check inputs for sanity
  int sumOfPowers = a_power.sum();
  if(sumOfPowers <= 0)
    {
      MayDay::Error("invalid power sum for monomialif");
    }
  for(int ipow = 0; ipow < SpaceDim-1; ipow++)
    {
      if(a_power[ipow] < 0)
        {
          MayDay::Error("invalid power for monomialif");
        }
    }

  m_inside = a_inside;
  m_power  = a_power; 
}
/*******/
Real 
MonomialIF::
value(const IndexTM<int,SpaceDim> & a_partialDerivative,
      const IndexTM<Real,SpaceDim>& a_point) const
{
  Real retval;
  int order = a_partialDerivative.sum();
  if (order == 0)
    {
      retval = value(a_point);
    }
  else
    {
      //partial_k(x^p) = pchoosek*(x^(k-p))
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_partialDerivative[idir] < 0)
            {
              MayDay::Error("invalid partial derivative");
            }
        }
      //check to see if there are x or y derivatives
      bool xyDerivatives = false;
      for(int idir = 0; idir < SpaceDim-1; idir++)
        {
          if(a_partialDerivative[idir] > 0)
            {
              xyDerivatives = true;
            }
        }
      //first the z part of this
      // if the power is too high, return zero
      // form of this function is -z + x^p1 y^2
      Real zval;
      if(xyDerivatives || (a_partialDerivative[SpaceDim-1] > 1))
        {
          zval = 0;
        }
      else if(a_partialDerivative[SpaceDim-1] == 1)
        {
          zval = -1; 
        }
      else if(a_partialDerivative[SpaceDim-1] == 0)
        {
          //should be covered by order==0
          zval = -a_point[SpaceDim-1];
        }
      else
        {
          MayDay::Error("missed a case");
        }
      //now for the xy val;
      Real xyval = 0;
      bool skipxy = (a_partialDerivative[SpaceDim-1] > 0);
      for(int idir = 0; idir < SpaceDim-1; idir++)
        {
          if(a_partialDerivative[idir] > m_power[idir])
            {
              skipxy = true;
            }
        }

      if(!skipxy)
        {
          //last one is spacedim = zval
          IndexTM<Real, SpaceDim-1> sdminonePoint;
          IndexTM<int , SpaceDim-1> sdminonePower;
          IndexTM<int , SpaceDim-1> sdminoneDeriv;
          for(int idir = 0; idir < SpaceDim-1; idir++)
            {  
              sdminonePoint[idir]= a_point[idir];
              sdminonePower[idir]= m_power[idir];
              sdminoneDeriv[idir]= a_partialDerivative[idir];
            }
  
          IndexTM<int, SpaceDim-1> pmink = sdminonePower-sdminoneDeriv;
          
          Real pminkfact = pfactorial(pmink);
          Real pfact     = pfactorial(sdminonePower);
          Real pchoosek = pfact/pminkfact;

          //already know that power >= deriv
          Real xtothepmink = power(sdminonePoint, pmink);
  
          xyval = pchoosek*xtothepmink;
        }

      retval = zval + xyval;
    }
  //  pout() << "retval2  = " << retval << endl;
  return retval;
}

Real 
MonomialIF::
value(const RealVect& a_point) const
{
  IndexTM<Real,GLOBALDIM> pt;

  //check does SpaceDim = GLOBALDIM
  CH_assert(GLOBALDIM == SpaceDim);
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      pt[idir] = a_point[idir];
    }
  //  pout() << "value 1  = " << value(pt) << endl;
  return value(pt);
}
Real 
MonomialIF::
xToTheQ(const IndexTM<Real,SpaceDim>& a_point) const
{
  IndexTM<Real, SpaceDim-1> sdminonePoint;
  IndexTM<int , SpaceDim-1> sdminonePower;
  for(int idir = 0; idir < SpaceDim-1; idir++)
    {  
      sdminonePoint[idir]= a_point[idir];
      sdminonePower[idir]= m_power[idir];
    }
  Real xyval = power(sdminonePoint, sdminonePower);
  return xyval;
}
Real 
MonomialIF::
value(const IndexTM<Real,SpaceDim>& a_point) const
{
  Real zval = -a_point[SpaceDim-1];
  //last one is spacedim = zval
  Real xyval = xToTheQ(a_point);

  Real retval = xyval + zval;
  if(!m_inside)
    retval = -retval;
  
  //  pout() << "retval  = " << retval << endl;
  return retval;
}

BaseIF* 
MonomialIF::
newImplicitFunction() const
{
  MonomialIF* MonomialPtr = new MonomialIF(m_power,m_inside);
  return static_cast<BaseIF*>(MonomialPtr);
}

///the stuff below is for geometryshop testing                                          
///return int x^p dV for the vof
IndMomSpaceDim 
MonomialIF::
getUncutVolumeSection(const RealVect& a_xlo, 
                      const RealVect& a_xhi, 
                      const RealVect& a_xbar,
                      const Real    & a_dx) const
{
  CH_assert(SpaceDim==2);
  IndMomSpaceDim retval;

  //  int Q = m_power[0];
  Real x0 = a_xlo[0];
  Real x1 = a_xhi[0];
  Real y0 = a_xlo[1];
  Real y1 = a_xhi[1];
  for(MomItSpaceDim momit; momit.ok(); ++momit)
    {
      IndexTM<int, SpaceDim> p = momit();
      Real xbit = (pow(x1, p[0]+1)-pow(x0, p[0]+1))/(p[0]+1);
      Real ybit = (pow(y1, p[1]+1)-pow(y0, p[1]+1))/(p[1]+1);
      Real thismom = xbit*ybit;
      retval[p]    = thismom;
    }

  //this integral was written for xbar = 0
  //need to shift to cell center
  shiftToXBar(retval, a_xbar);

  return retval;
}
///
void 
MonomialIF::
shiftToXBar(IndMomSpaceDim& a_moment,
            const RealVect& a_xbar) const
{
  RvSpaceDim shiftitm;
  //this integral was written for xbar = 0
  //need to shift to cell center
  RealVect shift = -a_xbar;
  EBArith::convertToITM(shiftitm,shift);
  a_moment.shift(shiftitm);
}
///////
IndMomSpaceDim 
MonomialIF::
getCutVolumeSection(const RealVect & a_xlo, 
                    const RealVect & a_xhi, 
                    const RealVect & a_xbar,
                    const Real     & a_dx) const
{
  
  CH_assert(SpaceDim==2);
  IndMomSpaceDim retval;
  
  int Q = m_power[0];
  Real x0 = a_xlo[0];
  Real x1 = a_xhi[0];
  //  Real y0 = a_xlo[1];
  Real y1 = a_xhi[1];
  for(MomItSpaceDim momit; momit.ok(); ++momit)
    {
      IndexTM<int, SpaceDim> p = momit();
      int px = p[0]; 
      int py = p[1]; 
      Real rpow1   = px + 1;
      Real rpow2   = px + Q*(py+1) + 1;

      Real factor1 = pow(y1, py+1.)/((py+1.)*rpow1);
      Real factor2 = 1./((py+1.)*rpow2);


      Real bit1    = factor1*(pow(x1,rpow1) - pow(x0, rpow1));
      Real bit2    = factor2*(pow(x1,rpow2) - pow(x0, rpow2));
      Real thismom = bit1 - bit2;
      retval[p]    = thismom;
    }


  //this integral was written for xbar = 0
  //need to shift to cell center
  shiftToXBar(retval, a_xbar);
  return retval;
}
///////
IndMomSpaceDim 
MonomialIF::
getExactVolumeMoments(const VolIndex    & a_vof,
                      const Real        & a_dx) const
{
  IndMomSpaceDim retval;
  CH_assert(SpaceDim == 2);
  //  MayDay::Error("3D moments not implemented");

  RealVect vectdx = a_dx*RealVect::Unit;
  RealVect offset = 0.5*a_dx*RealVect::Unit;
  RealVect locorn = EBArith::getVoFLocation(a_vof, vectdx, -offset);
  RealVect hicorn = EBArith::getVoFLocation(a_vof, vectdx,  offset);
  RealVect middle = 0.5*(locorn + hicorn);

  Real x0 = locorn[0]; 
  Real y0 = locorn[1]; 
  Real x1 = hicorn[0]; 
  Real y1 = hicorn[1]; 
  IndexTM<Real, SpaceDim> pt0, pt1;
  EBArith::convertToITM(pt0, locorn);
  EBArith::convertToITM(pt1, hicorn);


  Real yvallo = xToTheQ(pt0);
  Real yvalhi = xToTheQ(pt1);
  //we eliminate some cases because this is monotonic
  if(yvalhi <= y0)
    {
      //we are above the curve
      retval = getUncutVolumeSection(locorn, hicorn, middle, a_dx);
    }
  else if(yvallo >= y1)
    {
      //we are below the curve
      retval.setToZero();
    }
  else if((y0 < yvallo) && (yvallo < y1))
    {
      //curve cuts the left line
      //ok so we hit the curve
      //this is the x value where the curve intersects the higher line
      Real xhi = pow(y1, 1.0/m_power[0]);
      Real xf = Min(xhi, x1);
      RealVect rvhicorn = hicorn;
      rvhicorn[0] = xf;
      retval = getCutVolumeSection(locorn, rvhicorn, middle, a_dx);
    }
  else 
    {
      //curve cuts the bottom line
      //this is the x value where the curve intersects the lower line
      Real xloline = pow(y0, 1.0/m_power[0]);      
      //this is the x value where the curve intersects the upper line
      Real xhiline = pow(y1, 1.0/m_power[0]);      

      Real xinterlo = Max(xloline, x0);
      Real xinterhi = Min(xhiline, x1);
      IndMomSpaceDim regSection, irrSection;
      {
        RealVect rvxlo = locorn;
        RealVect rvxhi = hicorn;
        rvxhi[0] = xinterlo;
        regSection = getUncutVolumeSection(rvxlo, rvxhi , middle, a_dx);
      }
      {
        RealVect rvxlo = locorn;
        RealVect rvxhi = hicorn;
        rvxlo[0] = xinterlo;
        rvxhi[0] = xinterhi;
        irrSection = getCutVolumeSection(  rvxlo, rvxhi , middle, a_dx);
      }
      retval =  regSection;
      retval += irrSection;
    }
  ///cut moments and uncut moments are already shifted.
  return retval;
}
  
void 
MonomialIF::
getTloThi(bool& a_noEB,
          Real& a_tlo, 
          Real& a_thi,
          const VolIndex   & a_vof,
          const Real       & a_dx) const
{

  RealVect vectdx = a_dx*RealVect::Unit;
  RealVect offset = 0.5*a_dx*RealVect::Unit;
  RealVect locorn = EBArith::getVoFLocation(a_vof, vectdx, -offset);
  RealVect hicorn = EBArith::getVoFLocation(a_vof, vectdx,  offset);
  RealVect middle = 0.5*(locorn + hicorn);
  Real tlo, thi;
  //  Real x0 = locorn[0]; 
  Real y0 = locorn[1]; 
  Real x1 = hicorn[0]; 
  Real y1 = hicorn[1]; 
  IndexTM<Real, SpaceDim> pt0, pt1;
  EBArith::convertToITM(pt0, locorn);
  EBArith::convertToITM(pt1, hicorn);


  Real yvallo = xToTheQ(pt0);
  Real yvalhi = xToTheQ(pt1);
  a_noEB = false;
  //we eliminate some cases because this is monotonic
  if(yvalhi <= y0)
    {
      //we are above the curve
      a_noEB = true;
      return;
    }
  else if(yvallo >= y1)
    {
      //we are below the curve
      a_noEB = true;
      return;
    }
  else if((y0 < yvallo) && (yvallo < y1))
    {
      a_noEB = false;
      //curve cuts the left line
      //ok so we hit the curve
      //this is the x value where the curve intersects the higher line
      Real xhi = pow(y1, 1.0/m_power[0]);
      Real xf = Min(xhi, x1);
      tlo = locorn[0];
      thi = xf;
    }
  else 
    {
      a_noEB = false;
      //curve cuts the bottom line
      //this is the x value where the curve intersects the lower line
      //curve cuts the bottom line
      //this is the x value where the curve intersects the lower line
      Real xlo = pow(y0, 1.0/m_power[0]);      
      if(xlo > x1)
        {
          MayDay::Error("seems I have missed a case");
        }

      Real xhiline = pow(y1, 1.0/m_power[0]);      
      Real xhi = Min(xhiline, x1);
      tlo = xlo;
      thi = xhi;
    }
  a_tlo = tlo;
  a_thi = thi;
}
///return int_eb x^p dA for the eb
IndMomSpaceDim 
MonomialIF::
getExactEBMoments(const VolIndex   & a_vof,
                  const Real       & a_dx) const
{
  IndMomSpaceDim retval;
  IndMomSpaceDim tempdb;
  CH_assert(SpaceDim == 2);
  Real tlo, thi;
  bool noEB;
  getTloThi(noEB, tlo, thi, a_vof, a_dx);
  if(noEB)
    {
      retval.setToZero();
      return retval;
    }
  OneDIntegrator integrator;
  int num_bits = 1;
  for(int iorder = 0; iorder < 2*CH_EBIS_ORDER; iorder++)
    {
      num_bits *= 2;
    }
  
  for(MomItSpaceDim momit; momit.ok(); ++momit)
    {
      RefCountedPtr<BaseOneDFunc> func(new MonomialEBMomentFunc(m_power, momit()));
                                   
      Real integralVal = integrator.integral(func, tlo, thi, num_bits);
      retval[momit()] = integralVal;
      tempdb[momit()] = integralVal; 
    }

  //cell centered location
  RealVect xbar =  EBArith::getVoFLocation(a_vof, a_dx*RealVect::Unit, RealVect::Zero);

  //this integral was written for xbar = 0
  //need to shift to cell center
  shiftToXBar(retval, xbar);
  tempdb = retval;
  return retval;
}
                                          
                                          
///return int_eb x^p n_i dA for the eb
Real 
MonomialIF::
getEBNormMomX(const Real        & a_tlo,
              const Real        & a_thi,
              const IvSpaceDim  & a_p) const
{
  CH_assert(SpaceDim == 2);
  Real Q = m_power[0];
  Real Px    = a_p[0];
  Real Py    = a_p[1];
  Real tpow = Px + Q*Py + Q;
  Real rvhi = Q*pow(a_thi, tpow)/tpow;
  Real rvlo = Q*pow(a_tlo, tpow)/tpow;
  Real retval = rvhi - rvlo;
  return retval;
}

Real 
MonomialIF::
getEBNormMomY(const Real        & a_tlo,
              const Real        & a_thi,
              const IvSpaceDim  & a_p) const
{
  CH_assert(SpaceDim == 2);
  Real Q = m_power[0];
  Real Px    = a_p[0];
  Real Py    = a_p[1];
  Real tpow = Px + Q*Py + 1;
  Real rvhi = -pow(a_thi, tpow)/tpow;
  Real rvlo = -pow(a_tlo, tpow)/tpow;
  Real retval = rvhi - rvlo;
  return retval;
}
///
Real
MonomialIF::
getZerothDerivativeOfXNormal(const Real& a_xloc) const
{
  //(%i13) xnorm:(-Q)*(x^(Q-1))*((1 + Q*Q*(x^(2*Q-2)))^(-1/2));
  //                                     Q - 1
  //                                    x      Q
  //(%o13)                      - ---------------------
  //                                    2 Q - 2  2
  //                              sqrt(x        Q  + 1)
  Real retval;
  Real Q = m_power[0];
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);

  Real numer = -Q*pow(X, Q-1);
  Real denom = sqrt(argum);
  retval = numer/denom;

  return retval;
}
///
Real
MonomialIF::
getFirstDerivativeOfXNormal(const Real& a_xloc) const
{
  Real retval;
  const int  numTerms = 2;
  Real numer[numTerms];
  Real denom[numTerms];
  //(%i14) firstdx:diff(xnorm,x);
  //                 3 Q - 4  3                 Q - 2
  //                x        Q  (2 Q - 2)      x      (Q - 1) Q
  //(%o14)          ---------------------- - ---------------------
  //                    2 Q - 2  2     3/2         2 Q - 2  2
  //                2 (x        Q  + 1)      sqrt(x        Q  + 1)
  //
  Real Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      //                               2
      //                            8 X               2
      //(%o3)                   ------------- - --------------
      //                            2     3/2           2
      //                        (4 X  + 1)      sqrt(4 X  + 1)
      
      Real X = a_xloc;
      Real argum = 1 + 4*X*X;
      numer[0] =  8*X*X;
      numer[1] =  -2;
      denom[0] =   pow(argum, 3./2.);
      denom[1] =   pow(argum, 1./2.);

      retval = sumAllDivisions(numer, denom, numTerms);
    }
  else if(Q == 3)
    {
      //(%i20) firstdxqthree:diff(xnormqthree,X);
      //                                5
      //                            54 X             6 X
      //(%o20)                  ------------- - --------------
      //                            4     3/2           4
      //                        (9 X  + 1)      sqrt(9 X  + 1)

      Real X = a_xloc;
      Real argum = 1 + 9*pow(X,4.);
      numer[0] =  54*pow(X,5.);
      numer[1] =  -6*X;
      denom[0] =   pow(argum, 3./2.);
      denom[1] =   pow(argum, 1./2.);

      retval = sumAllDivisions(numer, denom, numTerms);
      
    }
  else
    {
      Real X = a_xloc;
      Real argum = 1 + Q*Q*pow(X, 2*Q-2);
      numer[0] =  (2.*Q-2.)*pow(Q, 3)*pow(X, 3*Q-4);
      numer[1] = -Q*(Q - 1.)*pow(X, Q-2);
      denom[0] = 2*pow(argum, 3./2.);
      denom[1] =   pow(argum, 1./2.);


      retval = sumAllDivisions(numer, denom, numTerms);
    }
  return retval;
}
///
Real
getSecondDerivXNormalGen(const Real& a_xloc, const int& a_Q)
{
  Real retval;
  //(%i15) secondx:diff(firstdx,x);                                                       
  //          Q - 3                      3 Q - 5  3                                      
  //         x      (Q - 2) (Q - 1) Q   x        Q  (2 Q - 2) (3 Q - 4)                  
  //(%o15) - ------------------------ + -------------------------------                  
  //                2 Q - 2  2                  2 Q - 2  2     3/2                       
  //          sqrt(x        Q  + 1)         2 (x        Q  + 1)                           
  //                        3 Q - 5          3                5 Q - 7  5          2       
  //                       x        (Q - 1) Q  (2 Q - 2)   3 x        Q  (2 Q - 2)        
  //                     + ----------------------------- - ------------------------       
  //                              2 Q - 2  2     3/2            2 Q - 2  2     5/2        
  //                          2 (x        Q  + 1)           4 (x        Q  + 1)          
  //                                                                                     
  const int  numTerms = 4;
  Real numer[numTerms];
  Real denom[numTerms];

  Real Q = a_Q;
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);
 
  numer[0] =    -(pow(X,  Q-3))*(Q-2)*(Q-1)*Q;
  numer[1] =     (pow(X,3*Q-5))*(pow(Q,3))*(2*Q-2)*(3*Q-4);
  numer[2] =     (pow(X,3*Q-5))*(pow(Q,3))*(2*Q-2)*(Q-1);
  numer[3] = -3.*(pow(X,5*Q-7))*(pow(Q,5))*(2*Q-2)*(2*Q-2);
  denom[0] =    (pow(argum, 1./2.));
  denom[1] = 2.*(pow(argum, 3./2.));
  denom[2] = 2.*(pow(argum, 3./2.));
  denom[3] = 4.*(pow(argum, 5./2.));

  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
///
Real
getSecondDerivXNormalQ2(const Real& a_xloc)
{
  Real retval;
  Real X = a_xloc;
  Real argum = 1 + 4*X*X;
  //(%i4) seconddxqtwo:diff(firstdxqtwo,X);
  //                                                 3
  //                             24 X            96 X
  //(%o4)                    ------------- - -------------
  //                             2     3/2       2     5/2
  //                         (4 X  + 1)      (4 X  + 1)
  const int  numTerms = 2;
  Real numer[numTerms];
  Real denom[numTerms];

 
  numer[0] =    24*X;
  numer[1] =    -96*X*X*X;
  denom[0] =    (pow(argum, 3./2.));
  denom[1] =    (pow(argum, 5./2.));

  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
///
Real
getSecondDerivXNormalQ3(const Real& a_xloc)
{
  Real retval;
  Real X = a_xloc;
  Real argum = 1 + 9*pow(X,4);
  //(%i21) seconddxqthree:diff(firstdxqthree,X);
  //                                          4                8
  //                       6             378 X           2916 X
  //(%o21)         - -------------- + ------------- - -------------
  //                         4            4     3/2       4     5/2
  //                 sqrt(9 X  + 1)   (9 X  + 1)      (9 X  + 1)

  const int  numTerms = 3;
  Real numer[numTerms];
  Real denom[numTerms];

 
  numer[0] =    -6.;
  numer[1] =    378*pow(X,4);
  numer[2] =    -2916*pow(X,8);
  denom[0] =    (pow(argum, 1./2.));
  denom[1] =    (pow(argum, 3./2.));
  denom[2] =    (pow(argum, 5./2.));

  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
////
Real
MonomialIF::
getSecondDerivativeOfXNormal(const Real& a_xloc) const
{
  Real retval;
  int Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      retval = getSecondDerivXNormalQ2(a_xloc);
    }
  else if(Q == 3)
    {
      retval = getSecondDerivXNormalQ3(a_xloc);
    }
  else
    {
      retval = getSecondDerivXNormalGen(a_xloc, Q);
    }
  return retval;
}
///
Real
getThirdDerivXNormalGen(const Real& a_xloc, const int& a_Q) 
{
  Real Q = a_Q;
  Real retval;
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);
  const int  numTerms = 8;
  Real numer[numTerms];
  Real denom[numTerms];
  //(%i16) thirdx:diff(secondx,x);                                                               
  //          Q - 4                                                                             
  //         x      (Q - 3) (Q - 2) (Q - 1) Q                                                   
  //(%o16) - --------------------------------                                                   
  //                    2 Q - 2  2                                                              
  //              sqrt(x        Q  + 1)         
  //                                                 

  numer[0] =  -pow(X,Q-4)*(Q-3)*(Q-2)*(Q-1)*Q;
  denom[0] =  (pow(argum, 1./2.));

  //    3 Q - 6  3                                                                               
  //   x        Q  (2 Q - 2) (3 Q - 5) (3 Q - 4)                                                 
  // + -----------------------------------------                                                 
  //                2 Q - 2  2     3/2                                                           
  //            2 (x        Q  + 1)                                                             

  numer[1] =    (pow(X, 3*Q-6 ))*pow(Q,3)*(2*Q-2)*(3*Q-5)*(3*Q-4);
  denom[1] =  2.*(pow(argum, 3./2.));

  //    3 Q - 6          3                                                                      
  //   x        (Q - 1) Q  (2 Q - 2) (3 Q - 5)                                                  
  // + ---------------------------------------                                                  
  //               2 Q - 2  2     3/2                                                            
  //           2 (x        Q  + 1)                                                               
  //

  numer[2] =      (pow(X, 3*Q-6))*(Q-1)*pow(Q,3)*(2*Q-2)*(3*Q-5);
  denom[2] =  2.*(pow(argum, 3./2.));

  //    3 Q - 6                  3            
  //   x        (Q - 2) (Q - 1) Q  (2 Q - 2)  
  // + -------------------------------------  
  //              2 Q - 2  2     3/2         
  //          2 (x        Q  + 1)            
  //  
  numer[3] =      (pow(X, 3*Q-6 ))*(Q-2)*(Q-1)*pow(Q,3)*(2*Q-2);
  denom[3] =  2.*(pow(argum, 3./2.));
  //
  //       5 Q - 8  5          2                         
  //    3 x        Q  (2 Q - 2)  (5 Q - 7)               
  //  - ----------------------------------               
  //              2 Q - 2  2     5/2                     
  //          4 (x        Q  + 1)                        
  //
  numer[4] =  -3.*(pow(X, 5*Q-8 ))*pow(Q,5)*(2*Q-2)*(2*Q-2)*(5*Q-7);
  denom[4] =  4.*(pow(argum, 5./2.));  
  //
  //      5 Q - 8  5          2          
  //   3 x        Q  (2 Q - 2)  (3 Q - 4)
  // - ----------------------------------
  //             2 Q - 2  2     5/2      
  //         4 (x        Q  + 1)         
  //
  numer[5] =  -3.*(pow(X, 5*Q-8 ))*pow(Q,5)*pow(2*Q-2,2)*(3*Q-4);
  denom[5] =  4.*(pow(argum, 5./2.));
  //  
  //      5 Q - 8          5          2                     
  //   3 x        (Q - 1) Q  (2 Q - 2)                        
  // - --------------------------------                     
  //            2 Q - 2  2     5/2                            
  //        4 (x        Q  + 1)                             
  //
  numer[6] =  -3.*(pow(X, 5*Q-8 ))*(  Q-1)*pow(Q,5)*pow(2*Q-2,2);
  denom[6] =  4.*(pow(argum, 5./2.));
  //       7 Q - 10  7          3                                                                
  //   15 x         Q  (2 Q - 2)                                                                
  // + --------------------------                                                               
  //         2 Q - 2  2     7/2                                                                  
  //     8 (x        Q  + 1)                                                                     
  //                                                                                             
  numer[7] =  15.*(pow(X, 7*Q-10))*pow(Q,7)*pow(2*Q-2,3);
  denom[7] =  8.*(pow(argum, 7./2.));

  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}

///
Real
getThirdDerivXNormalQ2(const Real& a_xloc)
{
  Real retval;
  Real X = a_xloc;
  const int  numTerms = 3;
  Real numer[numTerms];
  Real denom[numTerms];
  //(%i5) thirddxqtwo:diff(seconddxqtwo,X);
  //                                         2                4
  //                      24            576 X           1920 X
  //(%o5)            ------------- - ------------- + -------------
  //                     2     3/2       2     5/2       2     7/2
  //                 (4 X  + 1)      (4 X  + 1)      (4 X  + 1)

  Real argum = 4*X*X + 1;
  numer[0] = 24;
  numer[1] = -576*X*X;
  numer[2] = 1920*pow(X,4);
  denom[0] = pow(argum, 3./2.);
  denom[1] = pow(argum, 5./2.);
  denom[2] = pow(argum, 7./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}

///
Real
getThirdDerivXNormalQ3(const Real& a_xloc)
{
  Real retval;
  const int  numTerms = 3;
  Real numer[numTerms];
  Real denom[numTerms];
  Real X = a_xloc;
  Real argum = 9*pow(X,4) + 1;

  //(%i22) thirddxqthree:diff(seconddxqthree,X);
  //                          3               7               11
  //                    1620 X         43740 X        262440 X
  //(%o22)           ------------- - ------------- + -------------
  //                     4     3/2       4     5/2       4     7/2
  //                 (9 X  + 1)      (9 X  + 1)      (9 X  + 1)

  numer[0] =  1620*pow(X,3);
  numer[1] = -43740*pow(X,7);
  numer[2] = 262440*pow(X,11);
  denom[0] = pow(argum, 3./2.);
  denom[1] = pow(argum, 5./2.);
  denom[2] = pow(argum, 7./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}
///
Real
MonomialIF::
getThirdDerivativeOfXNormal(const Real& a_xloc) const
{
  Real retval;
  int Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      retval = getThirdDerivXNormalQ2(a_xloc);
    }
  else if(Q == 3)
    {
      retval = getThirdDerivXNormalQ3(a_xloc);
    }
  else
    {
      retval = getThirdDerivXNormalGen(a_xloc, Q);
    }
  return retval;
}
///
Real
getFourthDerivXNormalGen(const Real& a_xloc, const int& a_Q) 
{
  Real Q = a_Q;
  Real X = a_xloc;

  const int  numTerms = 16;
  Real numer[numTerms];
  Real denom[numTerms];
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);

  // (%i17) fourthx:diff(thirdx,x);                                                           
  //           Q - 5                                                                          
  //          x      (Q - 4) (Q - 3) (Q - 2) (Q - 1) Q                                        
  // (%o17) - ----------------------------------------                                        
  //                         2 Q - 2  2                                                       
  //                   sqrt(x        Q  + 1)                                                  
  numer[0] =     -(pow(X,  Q-5))*(Q-4)*(Q-3)*(Q-2)*(Q-1)*Q;
  denom[0] =       (pow(argum, 1./2.));

  //     3 Q - 7  3                                                                           
  //    x        Q  (2 Q - 2) (3 Q - 6) (3 Q - 5) (3 Q - 4)                                   
  //  + ---------------------------------------------------                                   
  //                      2 Q - 2  2     3/2                                                  
  //                  2 (x        Q  + 1)                                                     

  numer[1] =      (pow(X,3*Q-7))*pow(Q,3)*(2*Q-2)*(3*Q-6)*(3*Q-5)*(3*Q-4);
  denom[1] =    2.*(pow(argum, 3./2.));

  //     3 Q - 7          3                                                                   
  //    x        (Q - 1) Q  (2 Q - 2) (3 Q - 6) (3 Q - 5)                                     
  //  + -------------------------------------------------                                     
  //                     2 Q - 2  2     3/2                                                   
  //                 2 (x        Q  + 1)                                                      

  numer[2] =      (pow(X,3*Q-7))*(Q-1)*pow(Q,3)*(2*Q-2)*(3*Q-6)*(3*Q-5);
  denom[2] =    2.*(pow(argum, 3./2.));

  //     3 Q - 7                  3                                                           
  //    x        (Q - 2) (Q - 1) Q  (2 Q - 2) (3 Q - 6)                                       
  //  + -----------------------------------------------                                       
  //                    2 Q - 2  2     3/2                                                    
  //                2 (x        Q  + 1)                                                       

  numer[3] =     (pow(X,3*Q-7))*(Q-2)*(Q-1)*pow(Q,3)*(2*Q-2)*(3*Q-6);
  denom[3] = 2*pow(argum, 3./2.);

  //     3 Q - 7                          3                                                   
  //    x        (Q - 3) (Q - 2) (Q - 1) Q  (2 Q - 2)                                         
  //  + ---------------------------------------------                                         
  //                   2 Q - 2  2     3/2                                                     
  //               2 (x        Q  + 1)                                                        

  numer[4] =     (pow(X,3*Q-7))*(Q-3)*(Q-2)*(Q-1)*pow(Q,3)*(2*Q-2);
  denom[4] = 2*pow(argum, 3./2.);

  //       5 Q - 9  5          2                                                              
  //    3 x        Q  (2 Q - 2)  (5 Q - 8) (5 Q - 7)                                          
  //  - --------------------------------------------                                          
  //                   2 Q - 2  2     5/2                                                     
  //               4 (x        Q  + 1)                                                        

  numer[5] = -3.*(pow(X,5*Q-9))*pow(Q,5)*pow(2*Q-2,2)*(5*Q-8)*(5*Q-7);
  denom[5] = 4*pow(argum, 5./2.);
  //       5 Q - 9  5          2                                                              
  //    3 x        Q  (2 Q - 2)  (3 Q - 4) (5 Q - 8)                                          
  //  - --------------------------------------------                                          
  //                   2 Q - 2  2     5/2                                                     
  //               4 (x        Q  + 1)                                                        

  numer[6] = -3.*(pow(X,5*Q-9))*pow(Q,5)*pow(2*Q-2,2)*(3*Q-4)*(5*Q-8);
  denom[6] = 4*pow(argum, 5./2.);

  //       5 Q - 9          5          2                                                      
  //    3 x        (Q - 1) Q  (2 Q - 2)  (5 Q - 8)                                            
  //  - ------------------------------------------                                            
  //                  2 Q - 2  2     5/2                                                      
  //              4 (x        Q  + 1)                                                         

  numer[7] = -3.*(pow(X,5*Q-9))*(  Q-1)*pow(Q,5)*(2*Q-2)*(2*Q-2)*(5*Q-8);
  denom[7] = 4*pow(argum, 5./2.);

  //       5 Q - 9  5          2                                                              
  //    3 x        Q  (2 Q - 2)  (3 Q - 5) (3 Q - 4)                                          
  //  - --------------------------------------------                                          
  //                   2 Q - 2  2     5/2                                                     
  //               4 (x        Q  + 1)                                                        

  numer[8] = -3.*(pow(X,5*Q-9))*pow(Q,5)*pow(2*Q-2,2)*(3*Q-5)*(3*Q-4);
  denom[8] = 4*pow(argum, 5./2.);

  //       5 Q - 9          5          2                                                      
  //    3 x        (Q - 1) Q  (2 Q - 2)  (3 Q - 5)                                            
  //  - ------------------------------------------                                            
  //                  2 Q - 2  2     5/2                                                      
  //              4 (x        Q  + 1)                                    
  numer[9] = -3.*(pow(X,5*Q-9))*(Q-1)*pow(Q,5)*pow(2*Q-2,2)*(3*Q-5);
  denom[9] = 4*pow(argum, 5./2.);  

  //       5 Q - 9                  5          2                                              
  //    3 x        (Q - 2) (Q - 1) Q  (2 Q - 2)                                               
  //  - ----------------------------------------                                              
  //                 2 Q - 2  2     5/2                                                       
  //             4 (x        Q  + 1)                                                          

  numer[10]= -3.*(pow(X,5*Q-9))*(Q-2)*(Q-1)*pow(Q,5)*pow(2*Q-2,2);
  denom[10]= 4*pow(argum, 5./2.);


  //        7 Q - 11  7          3           
  //    15 x         Q  (2 Q - 2)  (7 Q - 10)
  //  + -------------------------------------
  //               2 Q - 2  2     7/2        
  //           8 (x        Q  + 1)           
  //
  //

  numer[11] =   15.*(pow(X,7*Q-11))*pow(Q,7)*pow(2*Q-2, 3)*(7*Q-10);
  denom[11] =   8*pow(argum, 7./2.);
 
  //
  //       7 Q - 11  7          3                    
  //   15 x         Q  (2 Q - 2)  (5 Q - 7)          
  // + ------------------------------------          
  //              2 Q - 2  2     7/2                 
  //          8 (x        Q  + 1)                    

  numer[12] =   15.*(pow(X,7*Q-11))*pow(Q,7)*pow(2*Q-2, 3)*(5*Q-7);
  denom[12] =   8*pow(argum, 7./2.);

  //
  //        7 Q - 11  7          3           
  //    15 x         Q  (2 Q - 2)  (3 Q - 4) 
  //  + ------------------------------------ 
  //               2 Q - 2  2     7/2        
  //           8 (x        Q  + 1)           
  //
  //
  numer[13] =   15.*(pow(X,7*Q-11))*pow(Q,7)*pow(2*Q-2, 3)*(3*Q-4);
  denom[13] =   8*pow(argum, 7./2.);

  //      7 Q - 11          7          3             
  //  15 x         (Q - 1) Q  (2 Q - 2)              
  //+ ----------------------------------             
  //            2 Q - 2  2     7/2                   
  //        8 (x        Q  + 1)                      
  //
  //

  numer[14]=    15.*(pow(X,7*Q-11))*(Q-1)*pow(Q,7)*pow(2*Q-2, 3);
  denom[14]=    8*pow(argum, 7./2.);

  //         9 Q - 13  9          4                                                           
  //    105 x         Q  (2 Q - 2)                                                            
  //  - ---------------------------                                                           
  //           2 Q - 2  2     9/2                                                             
  //      16 (x        Q  + 1)                                                                

  numer[15]=  -105.*(pow(X,9*Q-13))*pow(Q,9)*pow(2*Q-2, 4);
  denom[15]=   16*pow(argum, 9./2.);


  Real retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
///
Real
getFourthDerivXNormalQ2(const Real& a_xloc)
{
  Real retval;
  const int  numTerms = 3;
  Real numer[numTerms];
  Real denom[numTerms];
  Real X = a_xloc;
  Real argum = 4*X*X + 1;
  //  (%i6) fourthdxqtwo:diff(thirddxqtwo,X);
  //                                           3               5
  //                     1440 X         19200 X         53760 X
  //(%o6)           - ------------- + ------------- - -------------
  //                      2     5/2       2     7/2       2     9/2
  //                  (4 X  + 1)      (4 X  + 1)      (4 X  + 1)

  numer[0] = -1440*X;
  numer[1] = 19200*pow(X,3);
  numer[2] = -53760*pow(X,5);
  denom[0] = pow(argum, 5./2.);
  denom[1] = pow(argum, 7./2.);
  denom[2] = pow(argum, 9./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}
///
Real
getFourthDerivXNormalQ3(const Real& a_xloc)
{
  Real retval;
  const int  numTerms = 4;
  Real numer[numTerms];
  Real denom[numTerms];
  Real X = a_xloc;
  Real argum = 1 + 9*pow(X,4);
  //(%i23) fourthdxqthree:diff(thirddxqthree,X);
  //                  2                6               10              14
  //            4860 X         393660 X       6823440 X      33067440 X
  //(%o23)   ------------- - ------------- + ------------- - -------------
  //             4     3/2       4     5/2       4     7/2       4     9/2
  //         (9 X  + 1)      (9 X  + 1)      (9 X  + 1)      (9 X  + 1)

  numer[0] =      4860*pow(X,2);
  numer[1] =   -393660*pow(X,6);
  numer[2] =   6823440*pow(X,10);
  numer[3] = -33067440*pow(X,14);
  denom[0] = pow(argum, 3./2.);
  denom[1] = pow(argum, 5./2.);
  denom[2] = pow(argum, 7./2.);
  denom[3] = pow(argum, 9./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}
///
Real
MonomialIF::
getFourthDerivativeOfXNormal(const Real& a_xloc) const
{
  Real retval;
  int Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      retval = getFourthDerivXNormalQ2(a_xloc);
    }
  else if(Q == 3)
    {
      retval = getFourthDerivXNormalQ3(a_xloc);
    }
  else
    {
      retval = getFourthDerivXNormalGen(a_xloc, Q);
    }
  return retval;
}
///
Real
MonomialIF::
getZerothDerivativeOfYNormal(const Real& a_xloc) const
{
  // (%i18) ynorm:((1 + Q*Q*(x^(2*Q-2)))^(-1/2));                                                       
  //                                        1                                                           
  // (%o18)                       ---------------------                                                 
  //                                    2 Q - 2  2                                                      
  //                              sqrt(x        Q  + 1)                                                 
  //                                                                                                    
  Real retval = 0;
  Real Q = m_power[0];
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);
  Real numer = 1;
  Real denom = sqrt(argum);
  retval = numer/denom;

  return retval;
}
///
Real
MonomialIF::
getFirstDerivativeOfYNormal(const Real& a_xloc) const
{
  // (%i19) firstdy:diff(ynorm,x);                                                                      
  //                               2 Q - 3  2                                                           
  //                              x        Q  (2 Q - 2)                                                 
  // (%o19)                     - ----------------------                                                
  //                                  2 Q - 2  2     3/2                                                
  //                              2 (x        Q  + 1)                                                   

  Real retval = 0;
  Real Q = m_power[0];
  Real X = a_xloc;

  if(Q <= 1)
    {
      retval = 0;
    }
  else if (Q == 2)
    {
      //(%i9)  firstdyqtwo:diff(ynormqtwo,X);
      //                                       4 X
      //(%o9)                           - -------------
      //                                      2     3/2
      //                                  (4 X  + 1)
      Real numer =  -4*X;
      Real argum = 1 + 4*X*X;
      Real denom = pow(argum, 3./2.);
      retval = numer/denom;

    }
  else if (Q==3)
    {
      //firstdyqthree:diff(ynormqthree,X);
      //                                          3
      //                                      18 X
      //(%o25)                          - -------------
      //                                      4     3/2
      //                                  (9 X  + 1)
      Real numer =  -18*pow(X,3);
      Real argum = 1 + 9*pow(X,4);
      Real denom = pow(argum, 3./2.);
      retval = numer/denom;
    }
  else
    {

      //(%i8) firstderivgeny:diff(ynormgen,X);
      //                              2            2 Q - 3       
      //                             Q  (2 Q - 2) X              
      //(%o8)                      - ----------------------      
      //                                 2  2 Q - 2     3/2      
      //                             2 (Q  X        + 1)         
      //

      Real argum =  Q*Q*pow(X, 2*Q-2) + 1;
      Real numer =  -Q*Q*(2*Q-2)*(pow(X, 2*Q-3));
      Real denom = 2*pow(argum, 3./2.);
      retval = numer/denom;
    }
  return retval;
}
///
Real
getSecondDerivYNormalGen(const Real& a_xloc, const int& a_Q)
{
  
  Real retval = 0;
  Real Q = a_Q;
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2.*Q-2.);
  //(%i20) seconddy:diff(firstdy,x);                                                                   
  //             4 Q - 6  4          2    2 Q - 4  2                                                   
  //          3 x        Q  (2 Q - 2)    x        Q  (2 Q - 3) (2 Q - 2)                               
  //(%o20)    ------------------------ - -------------------------------                               
  //               2 Q - 2  2     5/2            2 Q - 2  2     3/2                                    
  //           4 (x        Q  + 1)           2 (x        Q  + 1)                                       
  //                                                                                                   
  const int  numTerms = 2;
  Real numer[numTerms];
  Real denom[numTerms];
  numer[0] = 3.*(pow(X,4.*Q-6.))*(pow(Q,4.))*pow(2.*Q-2., 2.);
  numer[1] =   -(pow(X,2.*Q-4.))*(pow(Q,2.))*(2.*Q-3.)*(2*Q-2);
  denom[0] =  4.*pow(argum, 5./2.);
  denom[1] =  2.*pow(argum, 3./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}

///
Real
getSecondDerivYNormalQ2(const Real& a_xloc)
{
  
  Real retval = 0;
  Real X = a_xloc;
  Real argum = 1 + 4*X*X;
  //(%i10) seconddyqtwo:diff(firstdyqtwo,X);
  //                                 2
  //                             48 X              4
  //(%o10)                   ------------- - -------------
  //                             2     5/2       2     3/2
  //                         (4 X  + 1)      (4 X  + 1)


  const int  numTerms = 2;
  Real numer[numTerms];
  Real denom[numTerms];
  numer[0] = 48*X*X;
  numer[1] = -4;
  denom[0] = pow(argum, 5./2.);
  denom[1] = pow(argum, 3./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}
///
Real
getSecondDerivYNormalQ3(const Real& a_xloc)
{
  
  Real retval = 0;
  Real X = a_xloc;
  Real argum = 1 + 9*pow(X,4);
  //(%i26) seconddyqthree:diff(firstdyqthree,X);
  //                                  6               2
  //                             972 X            54 X
  // (%o26)                   ------------- - -------------
  //                              4     5/2       4     3/2
  //                          (9 X  + 1)      (9 X  + 1)
  const int  numTerms = 2;
  Real numer[numTerms];
  Real denom[numTerms];
  numer[0] = 972*pow(X,6);
  numer[1] = -54*pow(X,2);
  denom[0] = pow(argum, 5./2.);
  denom[1] = pow(argum, 3./2.);
  retval = sumAllDivisions(numer, denom, numTerms);
  return retval;
}
///
Real
MonomialIF::
getSecondDerivativeOfYNormal(const Real& a_xloc) const
{
  Real retval;
  int Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      retval = getSecondDerivYNormalQ2(a_xloc);
    }
  else if(Q == 3)
    {
      retval = getSecondDerivYNormalQ3(a_xloc);
    }
  else
    {
      retval = getSecondDerivYNormalGen(a_xloc, Q);
    }
  return retval;
}
///
Real
getThirdDerivYNormalGen(const Real& a_xloc,const int& a_Q) 
{
  //(%i10) thirdderivgeny:diff(secondderivgeny,X);                                       
  //          2                                2 Q - 5                                   
  //         Q  (2 Q - 4) (2 Q - 3) (2 Q - 2) X                                          
  //(%o10) - -----------------------------------------                                   
  //                      2  2 Q - 2     3/2                                             
  //                  2 (Q  X        + 1)                                                
  //      4          2            4 Q - 7 
  //   3 Q  (2 Q - 2)  (4 Q - 6) X        
  // + ---------------------------------- 
  //             2  2 Q - 2     5/2       
  //         4 (Q  X        + 1)          
  //
  //
  //      4                    2  4 Q - 7           
  //   3 Q  (2 Q - 3) (2 Q - 2)  X                  
  // + ----------------------------------           
  //             2  2 Q - 2     5/2                 
  //         4 (Q  X        + 1)                    
  //
  //       6          3  6 Q - 9                                                         
  //   15 Q  (2 Q - 2)  X                                                                
  // - -------------------------                                                         
  //        2  2 Q - 2     7/2                                                           
  //    8 (Q  X        + 1)                                                              

  Real retval = 0;
  Real Q = a_Q;
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);
  const int  numTerms = 4;
  Real numer[numTerms];
  Real denom[numTerms];
  numer[0] =     -pow(Q,2)*(2*Q-4)*(2*Q-3)*(2*Q-2)*(pow(X, 2.*Q-5.));
  numer[1] =   3.*pow(Q,4)*pow(2*Q-2,2)*(4*Q-6)   *(pow(X, 4.*Q-7.));
  numer[2] =   3.*pow(Q,4)*pow(2*Q-2,2)*(2*Q-3)   *(pow(X, 4.*Q-7.));
  numer[3] =           -15.*pow(Q,6)*pow(2*Q-2, 3)*(pow(X, 6.*Q-9.));
  denom[0] =      2*pow(argum, 3./2.);
  denom[1] =      4*pow(argum, 5./2.);
  denom[2] =      4*pow(argum, 5./2.);
  denom[3] =      8*pow(argum, 7./2.);
  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}

///
Real
getThirdDerivYNormalQ2(const Real& a_xloc)
{
  //                                                                                                   
  //(%i17) thirddyqtwo:diff(seconddyqtwo,X);
  //                                                 3
  //                             144 X          960 X
  //(%o17)                   ------------- - -------------
  //                             2     5/2       2     7/2
  //                         (4 X  + 1)      (4 X  + 1)                                                                          
  Real retval = 0;                                     // 
  Real X = a_xloc;
  Real argum = 1 + 4*X*X;
  const int  numTerms = 2;
  Real numer[numTerms];
  Real denom[numTerms];
  numer[0] =  144*X;
  numer[1] = -960*pow(X,3);
  denom[0] = pow(argum, 5./2.);
  denom[1] = pow(argum, 7./2.);
  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
///
Real
getThirdDerivYNormalQ3(const Real& a_xloc)
{
  // (%i27) thirddyqthree:diff(seconddyqthree,X);
  //                                            5               9
  //                       108 X          8748 X         87480 X
  // (%o27)          - ------------- + ------------- - -------------
  //                       4     3/2       4     5/2       4     7/2
  //                   (9 X  + 1)      (9 X  + 1)      (9 X  + 1)

  Real retval = 0;
  Real X = a_xloc;
  Real argum = 1 + 9*pow(X,4);
  const int  numTerms = 3;
  Real numer[numTerms];
  Real denom[numTerms];
  numer[0] =   -108*X;
  numer[1] =   8748*pow(X,5);
  numer[2] = -87480*pow(X,9);
  denom[0] = pow(argum, 3./2.);
  denom[1] = pow(argum, 5./2.);
  denom[2] = pow(argum, 7./2.);
  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
///
Real
MonomialIF::
getThirdDerivativeOfYNormal(const Real& a_xloc) const
{
  Real retval;
  int Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      retval = getThirdDerivYNormalQ2(a_xloc);
    }
  else if(Q == 3)
    {
      retval = getThirdDerivYNormalQ3(a_xloc);
    }
  else
    {
      retval = getThirdDerivYNormalGen(a_xloc, Q);
    }
  return retval;
}
///
Real
getFourthDerivYNormalGen(const Real& a_xloc, const int& a_Q)
{
  Real retval = 0;
  Real Q = a_Q;
  Real X = a_xloc;
  Real argum = 1 + Q*Q*pow(X, 2*Q-2);
  const int  numTerms = 8;
  Real numer[numTerms];
  Real denom[numTerms];
  //                                                                                                    
  // (%i22) fourthdy:diff(thirddy,x);                                                                   
  //           2 Q - 6  2                                                                               
  //          x        Q  (2 Q - 5) (2 Q - 4) (2 Q - 3) (2 Q - 2)                                       
  // (%o22) - ---------------------------------------------------                                       
  //                            2 Q - 2  2     3/2                                                      
  //                        2 (x        Q  + 1)                                                         
  //       4 Q - 8  4          2                                                                        
  //    3 x        Q  (2 Q - 2)  (4 Q - 7) (4 Q - 6)                                                    
  //  + --------------------------------------------                                                                                                                                                      
  //                   2 Q - 2  2     5/2                                                               
  //               4 (x        Q  + 1)                                                                  
  //       4 Q - 8  4                    2                                                              
  //    3 x        Q  (2 Q - 3) (2 Q - 2)  (4 Q - 7)                                                    
  //  + --------------------------------------------                                                    
  //                   2 Q - 2  2     5/2                                                               
  //               4 (x        Q  + 1)                                                                  
  //       4 Q - 8  4                              2                                                    
  //    3 x        Q  (2 Q - 4) (2 Q - 3) (2 Q - 2)                                                     
  //  + --------------------------------------------                                                    
  //                   2 Q - 2  2     5/2                                                               
  //               4 (x        Q  + 1)                                                                                                                                                                    
  //        6 Q - 10  6          3                 6 Q - 10  6          3                               
  //    15 x         Q  (2 Q - 2)  (6 Q - 9)   15 x         Q  (2 Q - 2)  (4 Q - 6)                     
  //  - ------------------------------------ - ------------------------------------                     
  //               2 Q - 2  2     7/2                     2 Q - 2  2     7/2                            
  //           8 (x        Q  + 1)                    8 (x        Q  + 1)                               
  //        6 Q - 10  6                    3        8 Q - 12  8          4                              
  //    15 x         Q  (2 Q - 3) (2 Q - 2)    105 x         Q  (2 Q - 2)                               
  //  - ------------------------------------ + ---------------------------                              
  //               2 Q - 2  2     7/2                 2 Q - 2  2     9/2                                
  //           8 (x        Q  + 1)               16 (x        Q  + 1)                                   
  //                                                                                                    
  numer[0] =     -(pow(X,2*Q-6 ))*pow(Q,2)*(2*Q-5)*(2*Q-4)*(2*Q-3)*(2*Q-2);
  numer[1] =   3.*(pow(X,4*Q-8 ))*pow(Q,4)*pow(2*Q-2,2)*(4*Q-7)*(4*Q-6);
  numer[2] =   3.*(pow(X,4*Q-8 ))*pow(Q,4)*(2*Q-3)*pow(2*Q-2,2)*(4*Q-7);
  numer[3] =   3.*(pow(X,4*Q-8 ))*pow(Q,4)*(2*Q-4)*(2*Q-3)*pow(2*Q-2,2);
  numer[4] = -15.*(pow(X,6*Q-10))*pow(Q,6)*pow(2*Q-2,3)*(6*Q-9);
  numer[5] = -15.*(pow(X,6*Q-10))*pow(Q,6)*pow(2*Q-2,3)*(4*Q-6);
  numer[6] = -15.*(pow(X,6*Q-10))*pow(Q,6)*pow(2*Q-2,3)*(2*Q-3);
  numer[7] = 105.*(pow(X,8*Q-12))*pow(Q,8)*pow(2*Q-2,4); 
  denom[0] =  2.*pow(argum,3./2.);
  denom[1] =  4.*pow(argum,5./2.);
  denom[2] =  4.*pow(argum,5./2.);
  denom[3] =  4.*pow(argum,5./2.);
  denom[4] =  8.*pow(argum,7./2.);
  denom[5] =  8.*pow(argum,7./2.);
  denom[6] =  8.*pow(argum,7./2.);
  denom[7] = 16.*pow(argum,9./2.);
  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}

///
Real
getFourthDerivYNormalQ2(const Real& a_xloc)
{
  Real retval = 0;
  Real X = a_xloc;
  Real argum = 4*X*X + 1;
  const int  numTerms = 3;
  Real numer[numTerms];
  Real denom[numTerms];
  //(%i18) fourthddyqtwo:diff(thirddyqtwo,X);
  //                                          2               4
  //                      144           5760 X         26880 X
  //(%o18)           ------------- - ------------- + -------------
  //                     2     5/2       2     7/2       2     9/2
  //                 (4 X  + 1)      (4 X  + 1)      (4 X  + 1)

  numer[0] = 144;
  numer[1] = -5760*X*X;
  numer[2] = 26880*pow(X,4);
  denom[0] = pow(argum, 5./2.);
  denom[1] = pow(argum, 7./2.);
  denom[2] = pow(argum, 9./2.);
  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}

///
Real
getFourthDerivYNormalQ3(const Real& a_xloc)
{
  Real retval = 0;
  Real X = a_xloc;
  Real argum = 1 + 9*pow(X,4);
  const int  numTerms = 4;
  Real numer[numTerms];
  Real denom[numTerms];
  // (%i28) fourthdyqthree:diff(thirddyqthree,X);
  //                                    4                8               12
  //                108          49572 X        1574640 X      11022480 X
  // (%o28)  - ------------- + ------------- - ------------- + -------------
  //               4     3/2       4     5/2       4     7/2       4     9/2
  //           (9 X  + 1)      (9 X  + 1)      (9 X  + 1)      (9 X  + 1)

  numer[0] = -108;
  numer[1] =    49572*pow(X,4);
  numer[2] = -1574640*pow(X,8);
  numer[3] = 11022480*pow(X,12);
  denom[0] = pow(argum, 3./2.);
  denom[1] = pow(argum, 5./2.);
  denom[2] = pow(argum, 7./2.);
  denom[3] = pow(argum, 9./2.);
  retval = sumAllDivisions(numer, denom, numTerms);

  return retval;
}
///
Real
MonomialIF::
getFourthDerivativeOfYNormal(const Real& a_xloc) const
{
  Real retval;
  int Q = m_power[0];
  if(Q <= 1)
    {
      retval = 0;
    }
  else if(Q == 2)
    {
      retval = getFourthDerivYNormalQ2(a_xloc);
    }
  else if(Q == 3)
    {
      retval = getFourthDerivYNormalQ3(a_xloc);
    }
  else
    {
      retval = getFourthDerivYNormalGen(a_xloc, Q);
    }
  return retval;
}
///
IndMomSpaceDim 
MonomialIF::
getExactEBNormalPartialDerivs(const VolIndex   & a_vof,
                              const Real       & a_dx,
                              const int        & a_ni) const
{
  IndMomSpaceDim retval;
  CH_assert(SpaceDim == 2);
  RealVect loc = EBArith::getVoFLocation(a_vof, a_dx*RealVect::Unit, RealVect::Zero);
  //(all these things are only a function of x)
  Real X = loc[0];
  for(MomItSpaceDim momit; momit.ok(); ++momit)
    {
      int xpart = momit()[0];
      int ypart = momit()[1];
      if(ypart > 0) //no y derivs 
        {
          retval[momit()]  = 0;
        }
      else if(xpart > 4)
        {
          //I have not written the higher derivs yet
          retval[momit()]  = -7;
        }
      else if((a_ni == 0) && (xpart == 0))
        {
          retval[momit()] = getZerothDerivativeOfXNormal(X);
        }
      else if((a_ni == 0) && (xpart == 1))
        {
          retval[momit()] = getFirstDerivativeOfXNormal(X);
        }
      else if((a_ni == 0) && (xpart == 2))
        {
          retval[momit()] = getSecondDerivativeOfXNormal(X);
        }
      else if((a_ni == 0) && (xpart == 3))
        {
          retval[momit()] = getThirdDerivativeOfXNormal(X);
        }
      else if((a_ni == 0) && (xpart == 4))
        {
          retval[momit()] = getFourthDerivativeOfXNormal(X);
        }
      else if((a_ni == 1) && (xpart == 0))
        {
          retval[momit()] = getZerothDerivativeOfYNormal(X);
        }
      else if((a_ni == 1) && (xpart == 1))
        {
          retval[momit()] = getFirstDerivativeOfYNormal(X);
        }
      else if((a_ni == 1) && (xpart == 2))
        {
          retval[momit()] = getSecondDerivativeOfYNormal(X);
        }
      else if((a_ni == 1) && (xpart == 3))
        {
          retval[momit()] = getThirdDerivativeOfYNormal(X);
        }
      else if((a_ni == 1) && (xpart == 4))
        {
          retval[momit()] = getFourthDerivativeOfYNormal(X);
        }
      else
        {
          MayDay::Error("missed a case");
        }
    }
  //because outward facing;
  retval *= -1;
  return retval;
}
///
IndMomSpaceDim 
MonomialIF::
getExactEBNormalMoments(const VolIndex   & a_vof,
                        const Real       & a_dx,
                        const int        & a_ni) const
{
  IndMomSpaceDim retval;
  CH_assert(SpaceDim == 2);
  Real tlo, thi;
  bool noEB;
  getTloThi(noEB, tlo, thi, a_vof, a_dx);
  if(noEB)
    {
      retval.setToZero();
      return retval;
    }
  for(MomItSpaceDim momit; momit.ok(); ++momit)
    {
      Real integralVal;
      if(a_ni == 0)
        {
          integralVal =  getEBNormMomX(tlo, thi, momit());
        }
      else if(a_ni == 1)
        {
          integralVal =  getEBNormMomY(tlo, thi, momit());
        }
      else
        {
          MayDay::Error("bogus ni");
        }
      retval[momit()] = integralVal;
    }

  //cell centered location
  RealVect xbar =  EBArith::getVoFLocation(a_vof, a_dx*RealVect::Unit, RealVect::Zero);

  //this integral was written for xbar = 0
  //need to shift to cell center
  shiftToXBar(retval, xbar);

  return retval;
}
                                          
///return int x^p  dA for the face
IndMomSDMinOne 
MonomialIF::
getExactFaceMoments(const FaceIndex  & a_face,
                    const Real       & a_dx) const
{
  IndMomSDMinOne retval;
  retval.setToZero();
  MayDay::Error("not implemented");
  return retval;
}
#include "NamespaceFooter.H"
