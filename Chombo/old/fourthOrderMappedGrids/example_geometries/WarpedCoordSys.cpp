#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "WarpedCoordSys.H"

#include "UsingNamespace.H"

WarpedCoordSys::WarpedCoordSys()
{
  m_twoPi = 8.0*atan(1.0);
  m_isDefined = false;
}

WarpedCoordSys::~WarpedCoordSys()
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

void WarpedCoordSys::scale(const RealVect& a_scale)
{
   Real value = 1.0 / m_twoPi;
   RealVect bound( D_DECL( value, value, value ) );
   if ( a_scale > bound || a_scale < RealVect::Zero )
   {
      MayDay::Error("Bad scale in WarpedCoordSys::scale");
   }
   m_scale = a_scale;

   m_rootDir = SpaceDim;
   for (int d=0; d<SpaceDim; d++)
   {
      if (m_scale[d]>0)
      {
         m_rootDir = d;
         break;
      }
   }
}

void WarpedCoordSys::relative_tolerance(const Real a_rtol)
{
   if ( a_rtol < 0 )
   {
      MayDay::Error("Bad value in WarpedCoordSys::relative_tolerance");
   }
   m_rtol = a_rtol;
}

void WarpedCoordSys::absolute_tolerance(const Real a_atol)
{
   if ( a_atol < 0 )
   {
      MayDay::Error("Bad value in WarpedCoordSys::absolute_tolerance");
   }
   m_atol = a_atol;
}

void WarpedCoordSys::maximum_iterations(const int a_imax)
{
   if ( a_imax < 0 )
   {
      MayDay::Error("Bad value in WarpedCoordSys::maximum_iterations");
   }
   m_imax = a_imax;
}


void
WarpedCoordSys::define(const DisjointBoxLayout& a_grids,
                        const ProblemDomain& a_domain,
                        const RealVect& a_cellSpacing,
                        const IntVect& a_ghostVect)
{
  FourthOrderCoordSys::define(a_grids, a_domain,
                              a_cellSpacing, a_ghostVect);
  m_isDefined = true;
}



void
WarpedCoordSys::regrid(const DisjointBoxLayout& a_newGrids)
{
  // can't call regrid unless we've already "gridded"
  CH_assert(m_isDefined);

  // for now, just re-do define. May want to make this more
  // efficient in the future
  define(a_newGrids, m_domain, m_dx, m_ghostVect);
}

Real
WarpedCoordSys::g( const Real xi, const RealVect& x ) const
{
   const Real delta = (xi - x[m_rootDir]) / m_scale[m_rootDir];
   Real sine_product = 1.0;
   for (int d=0; d<SpaceDim; d++)
   {
      Real eta = x[d];
      if (d!=m_rootDir)
         eta += m_scale[d] * delta;
      else
         eta = xi;
      sine_product *= sin( m_twoPi * eta );
   }
   return x[m_rootDir] - m_scale[m_rootDir] * sine_product;
}

RealVect
WarpedCoordSys::realCoord(const RealVect& a_Xi) const
{
   Real sine_product = 1.0;
   for (int d=0; d<SpaceDim; d++)
      sine_product *= sin( m_twoPi * a_Xi[d] );
   RealVect xyzLoc;
   D_TERM(xyzLoc[0] = a_Xi[0] + m_scale[0] * sine_product;,
          xyzLoc[1] = a_Xi[1] + m_scale[1] * sine_product;,
          xyzLoc[2] = a_Xi[2] + m_scale[2] * sine_product;)
   return xyzLoc;
}

RealVect
WarpedCoordSys::mappedCoord(const RealVect& a_x) const
{
   RealVect mappedXi;
   for (int d=0; d<m_rootDir; d++)
      mappedXi[d] = a_x[d];

   if (m_rootDir<SpaceDim)
   {

      // find root by fixed point iteration
      Real root = a_x[m_rootDir];
      Real residual = g( root, a_x ) - root;
      Real bound = m_rtol * fabs(residual) + m_atol;
      int count = 0;
      while ( fabs(residual) > bound && count < m_imax )
      {
         root += residual;
         residual = g( root, a_x ) - root;
         count++;
      }
      if (count==m_imax)
      {
         MayDay::Error("Convergence failure in WarpedCoordSys::mappedCoord iteration!");
      }

      mappedXi[m_rootDir] = root;
      Real delta = (mappedXi[m_rootDir] - a_x[m_rootDir]) / m_scale[m_rootDir];
      for (int d=m_rootDir+1; d<SpaceDim; d++)
         mappedXi[d] = a_x[d] + m_scale[d] * delta;
   }

   return mappedXi;
}


Real
WarpedCoordSys::pointwiseJ(const RealVect& a_X) const
{
   return FourthOrderCoordSys::pointwiseJ( a_X );
}


void
WarpedCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                        LevelData<FluxBox>& a_F)
{
   FourthOrderCoordSys::mappedGridDivergence(a_divF, a_F);
}



Real
WarpedCoordSys::dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
   RealVect xi = mappedCoord( a_X );
   Real product = 1.0;
   for (int d=0; d<SpaceDim; d++)
   {
      Real arg = m_twoPi * xi[d];
      product *= ( d == a_dirXi ) ? cos( arg ) : sin( arg );
   }

   Real value = m_twoPi * m_scale[a_dirX] * product;
   if (a_dirX == a_dirXi)
   {
      value += 1.0;
   }

   return value;
}

// -- begin factory implementations ---------------------------

// default volInterval is all directions
WarpedCoordSysFactory::WarpedCoordSysFactory(const ProblemDomain& a_baseDomain,
                                             const Vector<int>& a_vectRefRatios,
                                             const RealVect& a_baseDx,
                                             const RealVect& a_scale,
                                             const Real a_rtol,
                                             const Real a_atol,
                                             const int a_imax)

  : m_volInterval(0,SpaceDim-1)
{
  int numLevels = a_vectRefRatios.size() +1;

  m_scale = a_scale;
  m_rtol = a_rtol;
  m_atol = a_atol;
  m_imax = a_imax;

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
WarpedCoordSysFactory::getCoordSys(const DisjointBoxLayout& a_grids,
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
      MayDay::Error("Invalid level in WarpedCoordSysFactory::getCoordSys");
    }

  WarpedCoordSys* newCSPtr = new WarpedCoordSys;

  newCSPtr->m_volInterval = m_volInterval;

  newCSPtr->scale(m_scale);
  newCSPtr->relative_tolerance(m_rtol);
  newCSPtr->absolute_tolerance(m_atol);
  newCSPtr->maximum_iterations(m_imax);

  newCSPtr->define(a_grids, m_vectDomains[level],
                   m_dxVect[level], a_ghostVect);


  return static_cast< CoordSys<FArrayBox,FluxBox>* >(newCSPtr);

}


