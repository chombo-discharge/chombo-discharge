#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianSphereIF.H"
#include "Factorial.H"
#include "NamespaceHeader.H"
///For these derivatives, let us be thankful for 
/// A) the existence of maxima
/// B) exponentials having uncomplicated derivatives
/// and most importantly
/// C) The equality of mixed derivs.
//// In all these cases:
/// f: C*exp(-b*(x*x + y*y + z*z));
///

bool isBogus(const Real& a_number)
{
  bool retval = false;
  if(a_number != a_number)
    retval = true;

  return retval;
}
Real 
GaussianSphereIF::
firstDeriv(const  IntVect& a_deriv,
           const RealVect& a_point) const
{
  CH_assert(a_deriv.sum() == 1);
  RealVect relLoc = a_point - m_center;
  Real rsq = relLoc.radSquared();

  //component of the point in the deriv direction;
  Real xd;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      CH_assert((a_deriv[idir]==0) || (a_deriv[idir]==1));
      if(a_deriv[idir] == 1)
        {
          xd = relLoc[idir];
        }
    }
  //  fx: diff(f,x)
  //                 2    2    2   
  //           - b (z  + y  + x )  
  // - 2 b x %e                   C

  Real retval = -2.*m_C*m_b*xd*exp(-m_b*rsq);
  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;
}
///
Real 
GaussianSphereIF::
secondDeriv(const  IntVect& a_deriv,
            const RealVect& a_point) const
{
  CH_assert(a_deriv.sum() == 2);
  RealVect relLoc = a_point - m_center;
  Real rsq = relLoc.radSquared();
  Real expval = exp(-m_b*rsq);
  //true for dxx dyy dzz, false for dxy dyz dxz
  //the two cases have different forms.
  bool doublex; 
  //component of the point in the deriv directions;
  Real xd[2];
  int iincr = 0;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      CH_assert((a_deriv[idir]>=0) && (a_deriv[idir]<=2));
      if(a_deriv[idir] == 1)
        {
          xd[iincr] = relLoc[idir];
          iincr++;
          doublex = false;
        }
      else if(a_deriv[idir] == 2)
        {
          xd[0] = relLoc[idir];
          xd[1] = relLoc[idir];
          doublex = true;
        }
    }
  Real retval = 0;
  if(doublex)
    {
      //(%i5) fxx: diff(fx,x);                                               
      //                          2    2    2                  2    2    2   
      //             2  2   - b (z  + y  + x )           - b (z  + y  + x )  
      //(%o5)     4 b  x  %e                   C - 2 b %e                   C
      //
      Real x = xd[0];
      retval = m_C*expval*(4.*m_b*m_b*x*x* - 2*m_b);
    }
  else
    {
      //(%i8) fxy: diff(fx,y);                                 
      //                                         2    2    2   
      //                           2       - b (z  + y  + x )  
      //(%o8)                   4 b  x y %e                   C
      //
      Real x = xd[0];
      Real y = xd[1];

      retval = m_C*expval*(4.*m_b*m_b*x*y);
    }
  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;

}
///
Real
GaussianSphereIF::
thirdDeriv(const  IntVect& a_deriv,
           const RealVect& a_point) const
{
  CH_assert(a_deriv.sum() == 3);
  RealVect relLoc = a_point - m_center;
  Real rsq = relLoc.radSquared();
  Real expval = exp(-m_b*rsq);
  //these are the different forms that the third deriv can take
  bool xxx = (a_deriv.max() == 3);
  bool xxy = (a_deriv.max() == 2);
#if CH_SPACEDIM==3
  bool xyz = (a_deriv.max() == 1);
#endif
  Real retval;
  if(xxx)
    {
      Real x = 0;
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 3)
            {
              found = true;
              x = relLoc[idir];
              break;
            }
        }
      if(!found) MayDay::Error("logic error gs3");
      ///(%i11) fxxx: diff(fxx,x);                                              
      ///                        2    2    2                      2    2    2   
      ///            2     - b (z  + y  + x )        3  3   - b (z  + y  + x )  
      ///(%o11)  12 b  x %e                   C - 8 b  x  %e                   C
      ///
      ///
      retval = m_C*expval*(12.*m_b*m_b*x - 8.*m_b*m_b*m_b*x*x*x);
    }
  else if(xxy)
    {
      Real x,y;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 2)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 1)
            {
              y = relLoc[idir];
              foundy = true;
            }
        }
      if((!foundx)||(!foundy)) MayDay::Error("logic error gs3.1");
      ///(%i16) fxxy: diff(fxx,y);                                              
      ///                      2    2    2                        2    2    2   
      ///          2     - b (z  + y  + x )        3  2     - b (z  + y  + x )  
      ///(%o16) 4 b  y %e                   C - 8 b  x  y %e                   C
      retval = m_C*expval*(4.*m_b*m_b*y - 8.*m_b*m_b*m_b*x*x*y);
    }
#if CH_SPACEDIM == 3
  else if(xyz)
    {
      CH_assert(a_deriv == IntVect::Unit);
      Real x = relLoc[0];
      Real y = relLoc[1];
      Real z = relLoc[2];
      ///(%i80) fxyz: diff(fxy,z);                                
      ///                                           2    2    2   
      ///                           3         - b (z  + y  + x )  
      ///(%o80)                - 8 b  x y z %e                   C
      retval = m_C*expval*(-8.*POW(m_b,3)*x*y*z);
    }
#endif
  else
    {
      MayDay::Error("missed a case gaussiansphere");
    }
  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;
}
///
Real
GaussianSphereIF::
fourthDeriv(const  IntVect& a_deriv,
            const RealVect& a_point) const
{
  CH_assert(a_deriv.sum() == 4);
  RealVect relLoc = a_point - m_center;
  Real rsq = relLoc.radSquared();
  Real expval = exp(-m_b*rsq);
  //these are the different forms that the fourth deriv can take
  bool xxxx = (a_deriv.max() == 4);
  bool xxxy = (a_deriv.max() == 3);
#if CH_SPACEDIM==2
  bool xxyy = (a_deriv.max() == 2);
#else 
  //there are more cases in 3D
  bool xxyy = ((a_deriv.max() == 2) && (a_deriv.min() == 0));
  bool xxyz = ((a_deriv.max() == 2) && (a_deriv.min() == 1));
#endif

  Real retval;
  if(xxxx)
    {
      Real x;
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 4)
            {
              x = relLoc[idir];
              found = true;
              break;
            }
        }
      if(!found) MayDay::Error("logic error gs4.0");
      ///(%i22) fxxxx: diff(fxxx,x);                                                    
      ///                        2    2    2                       2    2    2          
      ///           4  4   - b (z  + y  + x )         3  2   - b (z  + y  + x )         
      ///(%o22) 16 b  x  %e                   C - 48 b  x  %e                   C       
      ///                                                                 2    2    2   
      ///                                                       2   - b (z  + y  + x )  
      ///                                                 + 12 b  %e                   C

      retval = m_C*expval*(16.*POW(m_b, 4)*POW(x,4) - 48.*POW(m_b,3)*x*x + 12*m_b*m_b); 
    }
  else if(xxxy)
    {
      Real x,y;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 3)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 1)
            {
              y = relLoc[idir];
              foundy = true;
            }
        }
      if((!foundx)||(!foundy)) MayDay::Error("logic error gs4.05");
      ///(%i25) fxxxy: diff(fxxx,y);                                                
      ///                          2    2    2                        2    2    2   
      ///           4  3     - b (z  + y  + x )         3       - b (z  + y  + x )  
      ///(%o25) 16 b  x  y %e                   C - 24 b  x y %e                   C
      ///
      retval = m_C*expval*(16.*POW(m_b, 4)*POW(x,3)*y - 24.*POW(m_b,3)*x*y);
    }
  else if(xxyy)
    {
      Real xy[2];
      bool foundx[2];
      foundx[0] = false;
      foundx[1] = false;
      int iincr = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 2)
            {
              xy[iincr] = relLoc[idir];
              foundx[iincr] = true;
              iincr++;
            }
        }
      if((!foundx[0])||(!foundx[1])) MayDay::Error("logic error gs4.1");

      Real x = xy[0];
      Real y = xy[1];
      ///(%i32) fxxyy: diff(fxxy,y);                                                    
      ///                           2    2    2                      2    2    2        
      ///           4  2  2   - b (z  + y  + x )        3  2   - b (z  + y  + x )       
      ///(%o32) 16 b  x  y  %e                   C - 8 b  y  %e                   C     
      ///                                   2    2    2                   2    2    2   
      ///                      3  2   - b (z  + y  + x )        2   - b (z  + y  + x )  
      ///                 - 8 b  x  %e                   C + 4 b  %e                   C
      retval = m_C*expval*
        (16.*POW(m_b, 4)*x*x*y*y - 8.*POW(m_b,3)*y*y - 8.*POW(m_b,3)*x*x + 4.*m_b*m_b);

    }
#if CH_SPACEDIM==3
  else if(xxyz)
    {
      Real yz[2];
      Real x;
      bool foundx;
      foundx = false;
      int iincr = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 2)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 1)
            {
              yz[iincr] = relLoc[idir];
              iincr++;
            }
          else
            {
              MayDay::Error("logic error gs4.2");
            }
        }
      if(!foundx) MayDay::Error("logic error gs4.3");
      ///(%i84) fxxyz: diff(fxxy,z);                                                 
      ///                            2    2    2                       2    2    2   
      ///           4  2       - b (z  + y  + x )        3       - b (z  + y  + x )  
      ///(%o84) 16 b  x  y z %e                   C - 8 b  y z %e                   C
      ///
      Real y = yz[0];
      Real z = yz[1];
      
      retval = m_C*expval*(16.*POW(m_b, 4)*x*x*y*z - 8.*POW(m_b,3)*y*z);
    }
#endif
  else
    {
      MayDay::Error("missed a case in gaussiansphere 4th deriv");
    }
  
  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;
}
///
Real
GaussianSphereIF::
fifthDeriv(const  IntVect& a_deriv,
           const RealVect& a_point) const
{
  CH_assert(a_deriv.sum() == 5);
  RealVect relLoc = a_point - m_center;
  Real rsq = relLoc.radSquared();
  Real expval = exp(-m_b*rsq);
  //these are the different forms that the fifth deriv can take
  bool xxxxx = (a_deriv.max() == 5);
  bool xxxxy = (a_deriv.max() == 4);
#if CH_SPACEDIM==2
  bool xxxyy = (a_deriv.max() == 3);
#else 
  //there are more cases in 3D
  bool xxxyy = ((a_deriv.max() == 3) && (a_deriv.min() == 0));
  bool xxxyz = ((a_deriv.max() == 3) && (a_deriv.min() == 1));
  bool xxyyz =  (a_deriv.max() == 2);
#endif
  Real retval;
  if(xxxxx)
    {
      Real x;
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 5)
            {
              x = relLoc[idir];
              found = true;
              break;
            }
        }
      if(!found) MayDay::Error("logic error gs5.0");
      ///(%i35) fxxxxx: diff(fxxxx,x);                                                  
      ///                          2    2    2                        2    2    2       
      ///             5  5   - b (z  + y  + x )          4  3   - b (z  + y  + x )      
      ///(%o35) - 32 b  x  %e                   C + 160 b  x  %e                   C    
      ///                                                                 2    2    2   
      ///                                                     3     - b (z  + y  + x )  
      ///                                              - 120 b  x %e                   C
      retval = m_C*expval*
        (-32.*POW(m_b,5)*POW(x,5) + 160.*POW(m_b,4)*POW(x,3) -120.*POW(m_b,3)*x);
    }
  else if (xxxxy)
    {
      Real x,y;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 4)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 1)
            {
              y = relLoc[idir];
              foundy = true;
            }
        }
      if((!foundx)||(!foundy)) MayDay::Error("logic error gs5.05");
      ///(%i38) fxxxxy: diff(fxxxx,y);                                                  
      ///                            2    2    2                         2    2    2    
      ///             5  4     - b (z  + y  + x )         4  2     - b (z  + y  + x )   
      ///(%o38) - 32 b  x  y %e                   C + 96 b  x  y %e                   C 
      ///                                                                 2    2    2   
      ///                                                     3     - b (z  + y  + x )  
      ///                                               - 24 b  y %e                   C
      retval = m_C*expval*
        (-32.*POW(m_b, 4)*POW(x, 4)*y + 96.*POW(m_b,4)*x*x*y - 24.*POW(m_b,3)*y);

    }
  else if(xxxyy)
    {
      Real x,y;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 3)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 2)
            {
              y = relLoc[idir];
              foundy = true;
            }
        }
      if((!foundx)||(!foundy)) MayDay::Error("logic error gs5.1");

      ///(%i44) fxxxyy: diff(fxxxy,y);                                                  
      ///                             2    2    2                         2    2    2   
      ///             5  3  2   - b (z  + y  + x )         4    2   - b (z  + y  + x )  
      ///(%o44) - 32 b  x  y  %e                   C + 48 b  x y  %e                   C
      retval = m_C*expval*(-32.*POW(m_b, 5)*POW(x, 3)*y*y + 48.*POW(m_b,4)*x*y*y);
    }
#if CH_SPACEDIM==3
  else if(xxxyz)
    {
      Real yz[2];
      Real x;
      bool foundx = false;
      int iincr = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 3)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 1)
            {
              yz[iincr] = relLoc[idir];
              iincr++;
            }
          else
            {
              MayDay::Error("logic error gs5.2");
            }
        }
      if(!foundx) MayDay::Error("logic error gs5.3");
      Real y = yz[0];
      Real z = yz[1];
      ///(%i88) fxxxyz: diff(fxxxy,z);                                                  
      ///                           2    2    2                           2    2    2   
      ///           4         - b (z  + y  + x )         5  3       - b (z  + y  + x )  
      ///(%o88) 48 b  x y z %e                   C - 32 b  x  y z %e                   C

      retval = m_C*expval*(48.*POW(m_b, 4)*x*y*z - 32.*POW(m_b,5)*POW(x,3)*y*z);

    }
  else if(xxyyz)
    {
      Real xy[2];
      Real z;
      bool foundz = false;
      int iincr = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 2)
            {
              xy[iincr] = relLoc[idir];
              iincr++;
            }
          else if(a_deriv[idir] == 1)
            {
              z = relLoc[idir];
              foundz = true;
            }
          else
            {
              MayDay::Error("logic error gs5.4");
            }
        }
      if(!foundz) MayDay::Error("logic error gs5.5");
      Real x = xy[0];
      Real y = xy[1];
      ///(%i50) fxxyyz: diff(fxxyy,z);                                                  
      ///                               2    2    2                                     
      ///             5  2  2     - b (z  + y  + x )                                    
      ///(%o50) - 32 b  x  y  z %e                   C                                  
      ///                      2    2    2                         2    2    2          
      ///       4  2     - b (z  + y  + x )         4  2     - b (z  + y  + x )         
      /// + 16 b  y  z %e                   C + 16 b  x  z %e                   C       
      ///                  2    2    2                                                  
      ///      3     - b (z  + y  + x )                                                 
      /// - 8 b  z %e                   C                                               
      retval = m_C*expval*
        (- 32.*POW(m_b,5)*x*x*y*y*z 
         + 16.*POW(m_b,4)*y*y*z 
         + 16.*POW(m_b,4)*x*x*z
         -  8.*POW(m_b,3)*z);
    }
#endif
  else
    {
      MayDay::Error("missed a case in gaussiansphereif fifth deriv");
    }

  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;
}
///
Real
GaussianSphereIF::
sixthDeriv(const  IntVect& a_deriv,
           const RealVect& a_point) const
{
  CH_assert(a_deriv.sum() == 6);
  RealVect relLoc = a_point - m_center;
  Real rsq = a_point.radSquared();
  Real expval = exp(-m_b*rsq);
  //these are the different forms that the fourth deriv can take
  bool xxxxxx = (a_deriv.max() == 6);
  bool xxxxxy = (a_deriv.max() == 5);
#if CH_SPACEDIM==2
  bool xxxxyy = (a_deriv.max() == 4);
  bool xxxyyy = (a_deriv.max() == 3);
#else 
  //there are more cases in 3D
  bool xxxxyy = ((a_deriv.max() == 4) && (a_deriv.min() == 0));
  bool xxxxyz = ((a_deriv.max() == 4) && (a_deriv.min() == 1));
  bool xxxyyy = ((a_deriv.max() == 3) && (a_deriv.min() == 0));
  bool xxxyyz = ((a_deriv.max() == 3) && (a_deriv.min() == 1));
#endif
  Real retval;
  if(xxxxxx)
    {
      Real x;
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 6)
            {
              x = relLoc[idir];
              found = true;
              break;
            }
        }
      if(!found) MayDay::Error("logic error gs6.0");
      ///(%i53) fxxxxxx: diff(fxxxxx,x);                                                
      ///                        2    2    2                        2    2    2         
      ///           6  6   - b (z  + y  + x )          5  4   - b (z  + y  + x )        
      ///(%o53) 64 b  x  %e                   C - 480 b  x  %e                   C      
      ///                                 2    2    2                     2    2    2   
      ///                    4  2   - b (z  + y  + x )          3   - b (z  + y  + x )  
      ///             + 720 b  x  %e                   C - 120 b  %e                   C
      retval = m_C*expval*
        ( 64.*  POW(m_b,6)*POW(x,6) - 480.*POW(m_b,5)*POW(x,4)
          -720.*POW(m_b,4)*x*x      - 120.*POW(m_b,3));
      
    }
  else if(xxxxxy)
    {
      Real x,y;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 5)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 1)
            {
              y = relLoc[idir];
              foundy = true;
            }
        }
      if((!foundx)||(!foundy)) MayDay::Error("logic error gs6.05");
      ///(%i56) fxxxxxy: diff(fxxxxx,y);                                                
      ///                          2    2    2                          2    2    2     
      ///           6  5     - b (z  + y  + x )          5  3     - b (z  + y  + x )    
      ///(%o56) 64 b  x  y %e                   C - 320 b  x  y %e                   C  
      ///                                                                 2    2    2   
      ///                                                   4       - b (z  + y  + x )  
      ///                                            + 240 b  x y %e                   C
      retval = m_C*expval*
        ( 64.*  POW(m_b,6)*POW(x,5)*y - 320.*POW(m_b,5)*POW(x,3)*y
          +240.*POW(m_b,4)*x*y);

    }
  else if(xxxxyy)
    {
      Real x,y;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 4)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 2)
            {
              y = relLoc[idir];
              foundy = true;
            }
        }
      if((!foundx)||(!foundy)) MayDay::Error("logic error gs6.1");
      ///(%i44) fxxxyy: diff(fxxxy,y);                                                  
      ///                             2    2    2                         2    2    2   
      ///             5  3  2   - b (z  + y  + x )         4    2   - b (z  + y  + x )  
      ///(%o44) - 32 b  x  y  %e                   C + 48 b  x y  %e                   C
      ///                                2    2    2                      2    2    2   
      ///                   4  3   - b (z  + y  + x )         3     - b (z  + y  + x )  
      ///             + 16 b  x  %e                   C - 24 b  x %e                   C
      ///
      retval = m_C*expval*
        (-32.*POW(m_b,5)*POW(x,3)*y*y + 48.*POW(m_b,4)*x*y*y
         +16.*POW(m_b,4)*POW(x,3)     - 24.*POW(m_b,3)*x);

    }
  else if(xxxyyy)
    {
      Real xy[2];
      int numfound = 0;
      int iincr = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 3)
            {
              xy[iincr] = relLoc[idir];
              numfound++;
              iincr++;
            }
        }
      if(numfound != 2) MayDay::Error("logic error gs 6.2");
      Real x = xy[0];
      Real y = xy[1];
      ///(%i71) fxxxyyy: diff(fxxxyy,y);                                                
      ///                           2    2    2                         2    2    2     
      ///           6  3  3   - b (z  + y  + x )         5    3   - b (z  + y  + x )    
      ///(%o71) 64 b  x  y  %e                   C - 96 b  x y  %e                   C  
      ///                             2    2    2                         2    2    2   
      ///              5  3     - b (z  + y  + x )          4       - b (z  + y  + x )  
      ///        - 96 b  x  y %e                   C + 144 b  x y %e                   C
      retval = m_C*expval*
        (64. *POW(m_b,6)*POW(x,3)*POW(y,3) -  96.*POW(m_b,5)*x*POW(y,3)
         -96.*POW(m_b,5)*POW(x,3)*y        + 144.*POW(m_b,4)*x*y);
    }
#if CH_SPACEDIM==3
  else if(xxxyyz)
    {
      Real x,y,z;
      bool foundx=false, foundy=false, foundz= false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 3)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if(a_deriv[idir] == 2)
            {
              y = relLoc[idir];
              foundy = true;
            }
          else if(a_deriv[idir] == 1)
            {
              z = relLoc[idir];
              foundz = true;
            }
          else
            {
              MayDay::Error("logic error gs6.4");
            }
        }
      if((!foundx)||(!foundy)||(!foundz)) MayDay::Error("logic error gs6.5");
      ///(%i74) fxxxyyz: diff(fxxxyy,z);                                              
      ///                             2    2    2                                     
      ///           6  3  2     - b (z  + y  + x )                                    
      ///(%o74) 64 b  x  y  z %e                   C                                  
      ///                        2    2    2                         2    2    2      
      ///       5    2     - b (z  + y  + x )         5  3     - b (z  + y  + x )     
      /// - 96 b  x y  z %e                   C - 32 b  x  z %e                   C   
      ///                     2    2    2                                             
      ///       4       - b (z  + y  + x )                                            
      /// + 48 b  x z %e                   C                                          
      ///
      retval = m_C*expval*
        (64.  *POW(m_b,6)*POW(x,3)*y*y*z 
         -96. *POW(m_b,5)*x*y*y*z
         - 32.*POW(m_b,5)*POW(x,3)*z        
         + 48.*POW(m_b,4)*x*z);
    }
  else if(xxxxyz)
    {
      Real x,y,z;
      bool foundx=false, foundy=false, foundz= false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_deriv[idir] == 4)
            {
              x = relLoc[idir];
              foundx = true;
            }
          else if((a_deriv[idir] == 1) && !foundy)
            {
              y = relLoc[idir];
              foundy = true;
            }
          else if(a_deriv[idir] == 1)
            {
              z = relLoc[idir];
              foundz = true;
            }
          else
            {
              MayDay::Error("logic error gs6.6");
            }
        }
      if((!foundx)||(!foundy)||(!foundz)) MayDay::Error("logic error gs6.7");
      ///(%i68) fxxxxyz: diff(fxxxxy,z);                                                
      ///                            2    2    2                                        
      ///           6  4       - b (z  + y  + x )                                       
      ///(%o68) 64 b  x  y z %e                   C                                     
      ///                              2    2    2                        2    2    2   
      ///             5  2       - b (z  + y  + x )         4       - b (z  + y  + x )  
      ///      - 192 b  x  y z %e                   C + 48 b  y z %e                   C
      retval = m_C*expval*
        (64.  *POW(m_b,6)*POW(x,4)*y*z 
         -192.*POW(m_b,5)*x*x*y*z
         + 48.*POW(m_b,4)*y*z);
    }
#endif
  else
    {
      MayDay::Error("missed a case in gaussiansphereif sixth deriv");
    }

  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;
}
///
Real 
GaussianSphereIF::
derivative(const  IntVect& a_deriv,
           const RealVect& a_point) const
{
  Real retval = 0;
  int order = a_deriv.sum();
  if(order == 0)
    {
      retval = value(a_point);
    }
  else if(order == 1)
    {
      retval = firstDeriv(a_deriv, a_point);
    }
  else if(order == 2)
    {
      retval = secondDeriv(a_deriv, a_point);
    }
  else if(order == 3)
    {
      retval = thirdDeriv(a_deriv, a_point);
    }
  else if(order == 4)
    {
      retval = fourthDeriv(a_deriv, a_point);
    }
  else if(order == 5)
    {
      retval = fifthDeriv(a_deriv, a_point);
    }
  else if(order == 6)
    {
      retval = sixthDeriv(a_deriv, a_point);
    }
  else
    {
      //six is not enough? really?  
      retval = 0;
    }
  if(isBogus(retval))
    {
      MayDay::Error("bogosity detected");
    }
  return retval;
}
////
GaussianSphereIF::
GaussianSphereIF(const Real&     a_sigma,
                 const RealVect& a_center,
                 const bool&     a_inside)
{
  // Remember the parameters
  m_sigma  = a_sigma;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the sigma squared and the normalizing factor
  m_sigma2 = m_sigma * m_sigma;
  m_normal = 1.0 / (m_sigma * sqrt(2.0*M_PI));

  m_b = 1./(2.*m_sigma2);
  m_C = m_normal;
}
////
GaussianSphereIF::
GaussianSphereIF(const GaussianSphereIF& a_inputIF)
{
  // Remember the parameters
  m_sigma  = a_inputIF.m_sigma;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;

  m_sigma2 = a_inputIF.m_sigma2;
  m_normal = a_inputIF.m_normal;

  m_b      = a_inputIF.m_b;
  m_C      = a_inputIF.m_C;
}
////
GaussianSphereIF::~GaussianSphereIF()
{
}
////
void 
GaussianSphereIF::
GetParams(Real&     a_sigma,
          RealVect& a_center,
          bool&     a_inside) const
{
  // Copy parameter information over
  a_sigma  = m_sigma;
  a_center = m_center;
  a_inside = m_inside;
}
////
void 
GaussianSphereIF::
SetParams(const Real&     a_sigma,
          const RealVect& a_center,
          const bool&     a_inside)
{
  // Set parameter information
  m_sigma  = a_sigma;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the sigma squared and the normalizing factor
  m_sigma2 = m_sigma * m_sigma;
  m_normal = 1.0 / (m_sigma * sqrt(2.0*M_PI));

  m_b = 1./(2.*m_sigma2);
  m_C = m_normal;
}
////
Real 
GaussianSphereIF::
value(const RealVect& a_point) const
{
  Real retval;

  // The distance squared for m_center to a_point
  Real distance2;

  // Compute the distance squared
  distance2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real cur;
      cur = a_point[idir] - m_center[idir];

      distance2 += cur*cur;
    }

  retval = m_normal * exp(-distance2/(2*m_sigma2));

  // Change the sign to change inside to outside
  if (!m_inside)
    {
      retval = -retval;
    }

  return retval;
}
////
BaseIF* 
GaussianSphereIF::
newImplicitFunction() const
{
  GaussianSphereIF* spherePtr = new GaussianSphereIF(m_sigma,
                                                     m_center,
                                                     m_inside);

  return static_cast<BaseIF*>(spherePtr);
}
////
#include "NamespaceFooter.H"
