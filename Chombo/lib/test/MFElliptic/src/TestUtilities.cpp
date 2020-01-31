#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TestUtilities.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "SphereIF.H"
#include "PlaneIF.H"
#include "MFIndexSpace.H"
#include "DirichletPoissonDomainBC.H"
#include "LoadBalance.H"

Real ConstantDirichletBC::value (const RealVect &a_point,
                                 const RealVect &a_normal,
                                 const Real&     a_time,
                                 const int&      a_comp) const
{
  return m_constant;
}

Real TrilinearDirichletBC::value (const RealVect &a_point,
                                  const RealVect &a_normal,
                                  const Real&     a_time,
                                  const int&      a_comp) const
{
  Real retval = m_constant;
  for (int idir=0; idir<SpaceDim; ++idir)
    {
      retval += a_point[idir]*m_linearCoefficients[idir];
    }
  return retval;
}

Real TriquadraticDirichletBC::value (const RealVect &a_point,
                                     const RealVect &a_normal,
                                     const Real&     a_time,
                                     const int&      a_comp) const
{
  Real retval = m_constant;
  for (int idir=0; idir<SpaceDim; ++idir)
    {
      Real x = a_point[idir];
      retval += x*( m_linearCoefficients[idir] +
                    x*m_xxCoefficients[idir] );
      for (int jdir=idir+1; jdir<SpaceDim; ++jdir)
        {
          Real y = a_point[jdir];
          retval += x*y*m_xyCoefficients[jdir];
        }
    }
  return retval;
}

int readDomainInfo(int&       a_nlevels,
                   Box&       a_domain,
                   Real&      a_dx,
                   RealVect&  a_origin,
                   int&       a_maxBoxSize)
{
  // parse input file
  ParmParse pp;

  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);

  int maxLevel;
  pp.get("max_level", maxLevel);
  CH_assert(maxLevel >= 0);
  a_nlevels = maxLevel + 1;

  IntVect lo = IntVect::Zero;
  IntVect hi;

  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_origin[idir] = prob_lo[idir];
    }

  a_maxBoxSize = 32;
  pp.query("maxboxsize", a_maxBoxSize);

  return 0;
}

 void readSphereInfo(RealVect& a_center,
                     Real&     a_radius)
{
  // parse input file
  ParmParse pp;

  // get initial sphere parameters
  pp.get("radius", a_radius);

  Vector<Real> vectorCenter;
  pp.getarr("center", vectorCenter, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_center[idir] = vectorCenter[idir];
    }
} // end read from file

void readPlaneInfo(RealVect& a_normal,
                   RealVect& a_point)
{
  // parse input file
  ParmParse pp;

  Vector<Real> tempv;
  pp.getarr("point", tempv, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_point[idir] = tempv[idir];
    }

  pp.getarr("normal", tempv, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_normal[idir] = tempv[idir];
    }
}

void makeHierarchy(      Vector<DisjointBoxLayout>& a_dbl,
                   const ProblemDomain&             a_baseDomain,
                   const IntVectSet&                a_baseTags,
                   const Vector<int>&               a_refRatio,
                   const int                        a_nlevels,
                         int                        a_maxSize)
{
  a_dbl.resize(a_nlevels);
  Real fillRatio = 0.85;
  int blockingFactor= 4;
  int bufferSize =   2;
  BRMeshRefine regridder(a_baseDomain, a_refRatio, fillRatio, blockingFactor,
                         bufferSize, a_maxSize);

  Vector<Vector<Box> > oldGrids(a_nlevels, 1), newGrids(a_nlevels);
  oldGrids[0][0]=a_baseDomain.domainBox();
  for (int ilev=1; ilev<a_nlevels; ++ilev)
    {
      oldGrids[ilev][0] = coarsen(oldGrids[ilev-1][0], a_refRatio[ilev-1]);
    }

  regridder.regrid(newGrids, a_baseTags, 0, a_nlevels-1, oldGrids);

  Vector<int> procs;
  for (int i=0; i<a_nlevels; i++)
  {
    LoadBalance(procs, newGrids[i]);
    a_dbl[i] = DisjointBoxLayout(newGrids[i], procs);
  }
}

void zeroAdjacentToBoundary(LevelData<MFCellFAB>& a_data,
                            const ProblemDomain&  a_domain)
{
  DataIterator dit = a_data.dataIterator();
  const Box& d = a_domain.domainBox();


  for (dit.begin(); dit.ok(); ++dit)
    {
      MFCellFAB& mf = a_data[dit];
      for (int i=0 ;i<mf.numPhases(); ++i)
        {
          FArrayBox& regVal = mf.getPhase(i).getFArrayBox();
          for (int dir=0; dir<CH_SPACEDIM; ++dir)
            {
              if (!a_domain.isPeriodic(dir))
                {
                  Box adjHi = adjCellHi(d, dir, 1);  adjHi.shift(dir,-1);
                  Box adjLo = adjCellLo(d, dir, 1);  adjLo.shift(dir, 1);
                  regVal.setVal(0, adjHi & regVal.box(), 0, regVal.nComp());
                  regVal.setVal(0, adjLo & regVal.box(), 0, regVal.nComp());
                }
            }
        }
    }
}

void setValue(      LevelData<EBCellFAB>& a_phase,
              const BaseBCValue&          a_bc,
              const Box&                  a_domain,
              const RealVect&             a_dx,
              const RealVect&             a_origin,
                    bool                  a_useKappa)
{
  Real time = 0.;
  RealVect loc;
  IntVect  originIV = a_domain.smallEnd();
  DataIterator dit = a_phase.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      EBCellFAB& efab =  a_phase[dit];

      FArrayBox& fab  =  efab.getFArrayBox();
      ForAllX(Real, fab)
        {
          D_EXPR(loc[0]=iR+nR-nR+originIV[0], loc[1]=jR+originIV[1], loc[2]=kR+originIV[2]);
          loc += 0.5;
          loc *= a_dx;
          loc += a_origin;
          fabR = a_bc.value(loc, RealVect::Zero, time, 0);
        } EndFor ;
      if (a_useKappa)
        {
          const EBISBox& ebox = efab.getEBISBox();
          IntVectSet ivs = ebox.getIrregIVS(efab.getRegion());
          IVSIterator it(ivs);
          for (;it.ok(); ++it)
            {
              const IntVect iv = it();
              VolIndex vi(iv,0);
              if (ebox.bndryArea(vi) != 0)
                {
                  loc = iv;
                  loc += 0.5;
                  loc += ebox.centroid(vi);
                  loc *= a_dx;
                  loc += a_origin;
                  efab(vi, 0) = a_bc.value(loc, RealVect::Zero, time,0)*ebox.volFrac(vi);
                }
            }
        }
    }
}
