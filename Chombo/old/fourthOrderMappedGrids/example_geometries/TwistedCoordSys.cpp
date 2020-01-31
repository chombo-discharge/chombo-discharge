#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TwistedCoordSys.H"

#include "UsingNamespace.H"

TwistedCoordSys::TwistedCoordSys()
{
  m_Pi = 4.0*atan(1.0);
  m_isDefined = false;
}

TwistedCoordSys::~TwistedCoordSys()
{
  if (m_dNdXi.size() > 0)
    {
      for (int i=0; i<m_dNdXi.size(); i++)
        {
          if (m_dNdXi[i] != NULL)
            {
              delete m_dNdXi[i];
              m_dNdXi[i] = NULL;
            }
        }
    }
}

void
TwistedCoordSys::define(const DisjointBoxLayout& a_grids,
                        const ProblemDomain& a_domain,
                        const RealVect& a_cellSpacing,
                        const IntVect& a_ghostVect)
{
   FourthOrderCoordSys::define(a_grids, a_domain,
                               a_cellSpacing, a_ghostVect);
   m_isDefined = true;
}

void
TwistedCoordSys::regrid(const DisjointBoxLayout& a_newGrids)
{
  // can't call regrid unless we've already "gridded"
  CH_assert(m_isDefined);

  // for now, just re-do define. May want to make this more
  // efficient in the future
  define(a_newGrids, m_domain, m_dx, m_ghostVect);

}

// 14 Jun 2011:  Chombo/lib/src/BaseTools/Misc.H now has ipow
// inline Real ipow( const Real x, const int n )
// {
//    Real val( x );
//    for (int i=1; i<n; i++) val *= x;
//    return val;
// }

inline Real f( const Real r, const Real a )
{
   double val = (1.0 + cos( a * r ));
   return 0.5 * ipow(val,3);
}

#if 1
Real sinc( const Real x )
{
   Real val = 1.0;
   if ( fabs(x) > 1.0e-16 )
      val = sin(x) / x;
   else
   {
      static const Real A[] =
      {
        -1./6.,1./20., -1./42., 1./72., -1./110.
      };
      Real xsq = ipow(x,2);
      for (int i=4; i>=0; i--) val *= (1 - A[i] * xsq);
   }
   return val;
}

inline Real dfdronr( const Real r, const Real a )
{
   Real val = (1.0 + cos( a * r ));
   return (-1.5 * sinc( a * r ) * ipow(a*val,2) );
}
#else
inline Real dfdronr( const Real r, const Real a )
{
   Real val = (1.0 + cos( a * r ));
   return (-1.5 * a * sin( a * r ) * ipow(val,2) / r );
}
#endif

RealVect
TwistedCoordSys::realCoord(const RealVect& a_Xi) const
{
   RealVect xi(a_Xi);
   xi -= 0.5;
   Real r = sqrt( ipow(xi[0],2) + ipow(xi[1],2) );
   Real alpha = (r<m_R) ? f( r, m_scale) : 0;
   Real beta = alpha * m_theta;
   Real cosb = cos( beta );
   Real sinb = sin( beta );
   RealVect xyzLoc;
   D_TERM(xyzLoc[0] = cosb * xi[0] + sinb * xi[1];,
          xyzLoc[1] = cosb * xi[1] - sinb * xi[0];,
          xyzLoc[2] = xi[2];)
   xyzLoc += 0.5;
   return xyzLoc;
}

RealVect
TwistedCoordSys::mappedCoord(const RealVect& a_x) const
{
   RealVect x(a_x);
   x -= 0.5;
   Real r = sqrt( ipow(x[0],2) + ipow(x[1],2) );
   Real alpha = (r<m_R) ? f( r, m_scale) : 0;
   Real beta = alpha * m_theta;
   Real cosb = cos( beta );
   Real sinb = sin( beta );
   RealVect mappedXi;
   D_TERM(mappedXi[0] = cosb * x[0] - sinb * x[1];,
          mappedXi[1] = cosb * x[1] + sinb * x[0];,
          mappedXi[2] = x[2];)
   mappedXi += 0.5;
   return mappedXi;
}


Real
TwistedCoordSys::pointwiseJ(const RealVect& a_X) const
{
   return FourthOrderCoordSys::pointwiseJ( a_X );
}


void
TwistedCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                        LevelData<FluxBox>& a_F)
{
  FourthOrderCoordSys::mappedGridDivergence(a_divF, a_F);
}



Real
TwistedCoordSys::dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
   RealVect xi = mappedCoord( a_X );
   xi -= 0.5;
   Real r = sqrt( ipow(xi[0],2) + ipow(xi[1],2) );
   Real alpha = (r<m_R) ? f( r, m_scale) : 0;
   Real beta = alpha * m_theta;
   Real cosb = cos( beta );
   Real sinb = sin( beta );
   Real factor = (r<m_R) ? m_theta * dfdronr( r, m_scale ) : 0;

   Real value = 0.0;
   if (a_dirX == 0)
   {
      Real y = cosb * xi[1] - sinb * xi[0];
      if (a_dirXi == 0)
      {
         value =  cosb + xi[0] * y * factor;
      }
      else if (a_dirXi == 1)
      {
         value =  sinb + xi[1] * y * factor;
      }
      // else value stays 0
   }
   else if (a_dirX == 1)
   {
      Real x = cosb * xi[0] + sinb * xi[1];
      if (a_dirXi == 0)
      {
         value = -sinb - xi[0] * x * factor;
      }
      else if (a_dirXi == 1)
      {
         value =  cosb - xi[1] * x * factor;
      }
      // else value stays 0
   }
   else if (a_dirX == 2)
   {
      if (a_dirXi == 2)
      {
         value = 1.0;
      }
   }
   else
   {
      MayDay::Error("Bad dirX in TwistedCoordSys::dXdXi");
   }

   return value;
}

// -- begin factory implementations ---------------------------

// default volInterval is all directions
TwistedCoordSysFactory::TwistedCoordSysFactory(const ProblemDomain& a_baseDomain,
                                               const Vector<int>& a_vectRefRatios,
                                               const RealVect& a_baseDx,
                                               const Real& a_R,
                                               const Real& a_theta)
  : m_volInterval(0,SpaceDim-1)
{
  int numLevels = a_vectRefRatios.size() +1;

  m_R = a_R;
  m_theta = a_theta;

  // single-level case (if a single invalid refRatio is passed in)
  if (numLevels == 2 && a_vectRefRatios[0] <= 0)
    {
      numLevels = 1;
    }

  m_vectDomains.resize(numLevels);
  m_dxVect.resize(numLevels);
  m_vectRefRatios = a_vectRefRatios;

  m_vectDomains[0] = a_baseDomain;
  m_dxVect[0] = a_baseDx;

  for (int lev=1; lev<numLevels; lev++)
    {
      m_vectDomains[lev] = m_vectDomains[lev-1];
      m_vectDomains[lev].refine(m_vectRefRatios[lev-1]);
      m_dxVect[lev] = m_dxVect[lev-1];
      m_dxVect[lev] /= m_vectRefRatios[lev-1];
    }

}



CoordSys<FArrayBox, FluxBox>*
TwistedCoordSysFactory::getCoordSys(const DisjointBoxLayout& a_grids,
                                      const ProblemDomain& a_levelDomain,
                                      const IntVect& a_ghostVect) const
{
  // first identify level
  int level = -1;
  for (int lev=0; lev<m_vectDomains.size(); lev++)
    {
      if (m_vectDomains[lev] == a_levelDomain)
        {
          level = lev;
        }
    }

  if (level < 0)
    {
      pout() << "attempted to match bad domain: " << a_levelDomain
             << endl;
      MayDay::Error("Invalid level in TwistedCoordSysFactory::getCoordSys");
    }

  TwistedCoordSys* newCSPtr = new TwistedCoordSys;

  newCSPtr->m_volInterval = m_volInterval;

  newCSPtr->radius(m_R);
  newCSPtr->twist(m_theta);

  newCSPtr->define(a_grids, m_vectDomains[level],
                   m_dxVect[level], a_ghostVect);


  return static_cast< CoordSys<FArrayBox,FluxBox>* >(newCSPtr);

}


