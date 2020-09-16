#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BaseDomainBC.H"
#include "ViscousTensorOpF_F.H"
#include "EBViscousTensorOpF_F.H"
#include "BaseBCFuncEval.H"
#include "DirichletPoissonEBBC.H"
#include "ViscousTensorOp.H"
#include "NamespaceHeader.H"

//called by MAC projector enforceFaceVel
void
BaseDomainBC::
enforceFaceVel(LevelData<EBFluxFAB>&    a_velocity,
               const DisjointBoxLayout& a_grids,
               const EBISLayout&        a_ebisl,
               const ProblemDomain&     a_domain,
               const RealVect&          a_dx,
               const Real&              a_time,
               const RealVect&          a_origin)
{
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          DirichletPoissonEBBC::s_velComp = idir;
          if (!a_domain.isPeriodic(idir))
            {
              EBFluxFAB& fluxVel = a_velocity[dit()];
              CH_assert(fluxVel.nComp() == 1);
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Box sideBox = adjCellBox(grid, idir, sit(), 1);
                  int ishift = -sign(sit());
                  sideBox.shift(idir, ishift);
                  IntVectSet ivsIrreg(sideBox);
                  FaceStop::WhichFaces stopCrit;
                  if (sit() == Side::Lo)
                    {
                      stopCrit = FaceStop::LoBoundaryOnly;
                    }
                  else
                    {
                      stopCrit = FaceStop::HiBoundaryOnly;
                    }

                  FaceIterator faceit(ivsIrreg ,ebgraph, idir, stopCrit);
                  for (faceit.reset(); faceit.ok(); ++faceit)
                    {
                      Real boundaryVel = 0.0;
                      getFaceVel(boundaryVel, faceit(), fluxVel, a_origin, a_dx, idir, 0, a_time,sit());
                      fluxVel[idir](faceit(), 0) = boundaryVel;
                    }
                }
            }
        }
    }
  a_velocity.exchange();
}

void
BaseDomainBC::
enforceFaceVel(LevelData<EBFluxFAB>&    a_velocity,
               const DisjointBoxLayout& a_grids,
               const EBISLayout&        a_ebisl,
               const ProblemDomain&     a_domain,
               const RealVect&          a_dx,
               const Real&              a_time,
               const RealVect&          a_origin,
               const int&               a_comp)
{
  DirichletPoissonEBBC::s_velComp = a_comp;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (!a_domain.isPeriodic(idir))
            {
              EBFluxFAB& fluxVel = a_velocity[dit()];
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Box sideBox = adjCellBox(grid, idir, sit(), 1);
                  int ishift = -sign(sit());
                  sideBox.shift(idir, ishift);
                  IntVectSet ivsIrreg(sideBox);
                  FaceStop::WhichFaces stopCrit;
                  if (sit() == Side::Lo)
                    {
                      stopCrit = FaceStop::LoBoundaryOnly;
                    }
                  else
                    {
                      stopCrit = FaceStop::HiBoundaryOnly;
                    }

                  FaceIterator faceit(ivsIrreg ,ebgraph, idir, stopCrit);
                  for (faceit.reset(); faceit.ok(); ++faceit)
                    {
                      Real boundaryVel = 0.0;
                      getFaceVel(boundaryVel, faceit(), fluxVel, a_origin, a_dx, idir, 0, a_time,sit());
                      fluxVel[idir](faceit(), 0) = boundaryVel;
                    }
                }
            }
        }
    }
  a_velocity.exchange();
}


#include "NamespaceFooter.H"
