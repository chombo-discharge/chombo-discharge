#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DShapedCoordSys.H"

#include "UsingNamespace.H"

DShapedCoordSys::DShapedCoordSys()
{
  m_Pi = 4.0*atan(1.0);
  m_isDefined = false;
}

DShapedCoordSys::~DShapedCoordSys()
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
DShapedCoordSys::define(const DisjointBoxLayout& a_grids,
                        const double beta, const double kappa,
                        const RealVect& origin,
                        const ProblemDomain& a_domain,
                        const RealVect& a_cellSpacing,
                        const IntVect& a_ghostVect)
{
  m_beta = beta;
  m_kappa = kappa;
  m_origin = origin;

  FourthOrderCoordSys::define(a_grids, a_domain,
                              a_cellSpacing, a_ghostVect);


#if 0
  // not sure what I need this for
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

  m_dNdXi.resize(SpaceDim, NULL);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_dNdXi[dir] = new LevelData<FluxBox>(a_grids,
                                            SpaceDim, a_ghostVect);
    }
#endif

  m_isDefined = true;
}



void
DShapedCoordSys::regrid(const DisjointBoxLayout& a_newGrids)
{
  // can't call regrid unless we've already "gridded"
  CH_assert(m_isDefined);

  // for now, just re-do define. May want to make this more
  // efficient in the future
  define(a_newGrids, m_beta, m_kappa, m_origin, m_domain, m_dx, m_ghostVect);

}

RealVect
DShapedCoordSys::realCoord(const RealVect& a_X) const
{

  RealVect xyzLoc;
  D_TERM(xyzLoc[0] = m_origin[0] + a_X[0]*cos(a_X[1]+m_beta*sin(a_X[1]));,
         xyzLoc[1] = m_origin[1] + a_X[0]*m_kappa*sin(a_X[1]);,
         xyzLoc[2] = a_X[2];)

    return xyzLoc;
}

RealVect
DShapedCoordSys::mappedCoord(const RealVect& a_x) const
{
  double x = a_x[0] - m_origin[0];
  double y = a_x[1] - m_origin[1];
  double r;
  double theta;
  double two_pi = 2.0*m_Pi;

  /*
    We need to solve for theta, since it cannot be expressed
    analytically in terms of x and y for nonzero m_beta.

    We use a Newton iteration.  The initial guess for theta is exact
    when delta = m_beta = 0, and it should also be pretty good for
    nonzero values of the latter.  The iteration will be performed
    in the domain -pi < theta <= pi.  Before returning, theta is
    converted to the interval 0 <= theta < 2pi.
  */
  theta = atan2(y,x);

  double st, ct, t1, st1, ct1;

  while (true)
  {
    st  = sin(theta);
    ct  = cos(theta);
    t1  = theta + m_beta*st;
    st1 = sin(t1);
    ct1 = cos(t1);

    double g = m_kappa*x*st - y*ct1;

    if (fabs(g) > 1.e-11)
    {
       double gprime = m_kappa*x*ct + y*st1*(1. + m_beta*ct);
       theta -= g/gprime;
    }
    else break;
  }

  if ( fabs(st) > 1.e-11)
  {
    r = y / (m_kappa * st);
  }
  else
  {
    r = x / ct1;
  }

  // Now shift theta to the range 0 <= theta < 2pi.
  if (y < 0) theta += two_pi;

  RealVect mappedXi;
  D_TERM(mappedXi[0] = r;,
         mappedXi[1] = theta;,
         mappedXi[2] = a_x[2];)

  return m_stretch*mappedXi;

}


void
DShapedCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                        LevelData<FluxBox>& a_F)
{
  FourthOrderCoordSys::mappedGridDivergence(a_divF, a_F);
}



Real
DShapedCoordSys::dXdXi(const RealVect& a_X, int a_dirX, int a_dirXi) const
{
  Real value = 0.0;
  // this is for test purposes

  RealVect mapped_coord = mappedCoord(a_X);

  double x     = a_X[0];
  double y     = a_X[1];
  double r     = mapped_coord[0];
  double theta = mapped_coord[1];

  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          // dx/dr = cos(theta + beta*sin(theta))
          value = cos(theta + m_beta*sin(theta));
        }
      else if (a_dirXi == 1)
        {
          // dx/dtheta = -rsin(theta + beta*sin(theta))*(1 + beta*cos(theta))
          value = -r*sin(theta + m_beta*sin(theta))*(1. + m_beta*cos(theta));
        }
      // else value stays 0
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          // dy/dr = kappa*sin(theta)
          value = m_kappa * sin(theta);
        }
      else if (a_dirXi == 1)
        {
          // dy/dtheta = r kappa cos(theta)
          value = r * m_kappa * cos(theta);
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
      MayDay::Error("Bad dirX in RZThetaCoordSys::dXdXi");
    }

  return value;
}



// -- begin factory implementations ---------------------------

// default volInterval is all directions
DShapedCoordSysFactory::DShapedCoordSysFactory(const ProblemDomain& a_baseDomain,
                                               const double delta, const double kappa,
                                               const RealVect& origin,
                                               const Vector<int>& a_vectRefRatios,
                                               const RealVect& a_baseDx,
                                               const RealVect& a_stretch)
  : m_volInterval(0,SpaceDim-1),
    m_beta(asin(delta)),
    m_kappa(kappa),
    m_origin(origin)
{
  int numLevels = a_vectRefRatios.size() +1;

  m_stretch = a_stretch;

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
DShapedCoordSysFactory::getCoordSys(const DisjointBoxLayout& a_grids,
                                      const ProblemDomain& a_levelDomain,
                                      const IntVect& a_ghostVect)
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
      MayDay::Error("Invalid level in DShapedCoordSysFactory::getCoordSys");
    }

  DShapedCoordSys* newCSPtr = new DShapedCoordSys;

  newCSPtr->m_volInterval = m_volInterval;
  newCSPtr->stretch(m_stretch);

  newCSPtr->define(a_grids, m_beta, m_kappa, m_origin, m_vectDomains[level],
                   m_dxVect[level], a_ghostVect);


  return static_cast< CoordSys<FArrayBox,FluxBox>* >(newCSPtr);

}


