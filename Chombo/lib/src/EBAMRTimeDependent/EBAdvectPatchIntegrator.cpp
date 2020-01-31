#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBAdvectPatchIntegrator.H"
#include "EBPatchGodunovF_F.H"
#include "EBLoHiCenter.H"
#include "EBArith.H"
#include "PolyGeom.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "CH_Timer.H"
#include "EBPatchAdvectF_F.H"
#include "NamespaceHeader.H"
int  EBAdvectPatchIntegrator::s_doingVel  = -1;
int  EBAdvectPatchIntegrator::s_curComp   = -1;

/******/
EBAdvectPatchIntegrator::
EBAdvectPatchIntegrator(const Box&            a_validBox,
                        const EBISBox&        a_ebisBox,
                        const IntVectSet&     a_coarseFineIVS,
                        const ProblemDomain&  a_domain,
                        const RealVect&       a_dx,
                        bool                  a_useSlopeLimiting)
{
  CH_TIME("EBAdvectPatchIntegrator::setValidBox");
  m_isBCSet = false;
  m_useLimiting = a_useSlopeLimiting;
  m_validBox = a_validBox;
  m_ebisBox  = a_ebisBox;
  m_domain   = a_domain;
  m_dx       = a_dx;
  m_cfivs    = a_coarseFineIVS;
  m_isMaxMinSet = false;
  //I know that a lot of these objects are slightly larger than they
  //need to be.   Be warned, however, that it is
  //really important that all these things are the same size.
  //This allows me to make all kinds of wacky stenciling assumptions down
  //the road.  Mess with these definitions at thy peril.
  m_validBoxG4  = grow(a_validBox, 4);
  m_validBoxG4 &= m_domain;
  IntVectSet        ivsIrregG4 =  m_ebisBox.getIrregIVS(m_validBoxG4);
  VoFIterator vofit(ivsIrregG4,   m_ebisBox.getEBGraph());
  m_irregVoFs = vofit.getVector();

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVectSet irregIVSPlus, irregIVSMinu;
      EBArith::computeCoveredFaces(m_coveredFacePlusG4[idir],
                                   m_coveredSetsPlusG4[idir],
                                   irregIVSPlus,
                                   idir, Side::Hi, 
                                   m_ebisBox, m_validBoxG4);
      EBArith::computeCoveredFaces(m_coveredFaceMinuG4[idir],
                                   m_coveredSetsMinuG4[idir],
                                   irregIVSMinu,
                                   idir, Side::Lo, 
                                   m_ebisBox, m_validBoxG4);
    }
}
/***/
void
EBAdvectPatchIntegrator::
riemann(EBFaceFAB&       a_primGdnv,
        const EBCellFAB& a_primLeft,
        const EBCellFAB& a_primRigh,
        const int&       a_faceDir,
        const Box&       a_box)
{
  CH_TIME("EBADvectPatchIntegrator::riemann");
  //these two need to be set
  CH_assert((s_doingVel == 1) || (s_doingVel == 0));
  //  CH_assert((s_curComp >= 0) && (s_curComp < SpaceDim);
  CH_assert(s_curComp >= 0);
  CH_assert(m_isVelSet);

  // const EBCellFAB& normalVel = *m_normalVelPtr;
  const EBFluxFAB& advectVel = *m_advectionVelPtr;
  const EBFaceFAB& advectVelFaceFAB = advectVel[a_faceDir];
  int nPrim = 1;

  Box cellBox = enclosedCells(a_box);

  //explicitly cast the left and right states to modifiable references
  //because we need to shift them to faces and then shift them back.
  BaseFab<Real>& regPrimRigh = (BaseFab<Real>&)a_primRigh.getSingleValuedFAB();
  BaseFab<Real>& regPrimLeft = (BaseFab<Real>&)a_primLeft.getSingleValuedFAB();
  BaseFab<Real>& regPrimGdnv = a_primGdnv.getSingleValuedFAB();
  // const BaseFab<Real>& regNormalVel = normalVel.getSingleValuedFAB();
  const BaseFab<Real>& regAdvectVel = advectVelFaceFAB.getSingleValuedFAB();

  //find the regular part of qgdnv.
  FORT_ADVECTRIEMANN(CHF_BOX(a_box),
                     CHF_FRA(regPrimGdnv),
                     CHF_CONST_FRA(regPrimLeft),
                     CHF_CONST_FRA(regPrimRigh),
                     // CHF_CONST_FRA(regNormalVel),
                     CHF_CONST_FRA1(regAdvectVel, 0),//only 1 (normal) component at face
                     CHF_CONST_INT(a_faceDir),
                     CHF_CONST_INT(nPrim),
                     CHF_CONST_INT(s_curComp),
                     CHF_CONST_INT(s_doingVel));

  //the box sent into this is face-centered.
  //we need to use the cell-centered one it surrounds.
  Vector<FaceIndex> multiFaces =  m_ebisBox.getEBGraph().getMultiValuedFaces(a_faceDir, cellBox);
  for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (!face.isBoundary())
        {
          VolIndex vofl = face.getVoF(Side::Lo);
          VolIndex vofr = face.getVoF(Side::Hi);

          Real velFace = advectVelFaceFAB(face, 0);
          for (int ivar = 0; ivar < nPrim; ivar++)
            {
              if (velFace > EBAPI_TOL)
                {
                  a_primGdnv(face, ivar) = a_primLeft(vofl, ivar);
                }
              else if (velFace < -EBAPI_TOL)
                {
                  a_primGdnv(face, ivar) = a_primRigh(vofr, ivar);
                }
              else
                {
                  if ((s_doingVel == 1) && (s_curComp == a_faceDir))
                    {
                      a_primGdnv(face, ivar) = 0.0;
                    }
                  else
                    {
                      a_primGdnv(face, ivar) = 0.5*(a_primRigh(vofr, ivar) + a_primLeft(vofl, ivar));
                    }
                }
            }

        }
    }
}
/*****************************/
void
EBAdvectPatchIntegrator::
riemann(BaseIVFAB<Real>&        a_coveredPrim,
        const BaseIVFAB<Real>&  a_exteState,
        const EBCellFAB&        a_primState,
        const Vector<VolIndex>& a_vofset,
        const int&              a_faceDir,
        const Side::LoHiSide&   a_sd,
        const Box&       a_box)
{
  CH_TIME("EBADvectPatchIntegrator::riemannIrr");
  //this holds velocity at time n
  const EBCellFAB& normalVel = *m_normalVelPtr;
  int nPrim = 1;
  for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real  vel = normalVel(vof, a_faceDir);

          for (int ivar = 0; ivar < nPrim; ivar++)
            {
              Real leftState, righState;
              if (a_sd == Side::Lo)
                {
                  leftState = a_exteState(vof, ivar);
                  righState = a_primState(vof, ivar);
                }
              else
                {
                  leftState = a_primState(vof, ivar);
                  righState = a_exteState(vof, ivar);
                }

              if (vel > EBAPI_TOL)
                {
                  a_coveredPrim(vof, ivar) = leftState;
                }
              else if (vel < -EBAPI_TOL)
                {
                  a_coveredPrim(vof, ivar) = righState;
                }
              else
                {
                  a_coveredPrim(vof, ivar) = 0.5*(righState+leftState);
                }
            }
        }
    }
}
/******************/
void
EBAdvectPatchIntegrator::
consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const BaseIVFAB<Real>& a_ebIrregFlux,
                        const IntVectSet&      a_ivs)
{
  CH_TIME("EBADvectPatchIntegrator::consUndividedDivergence");
  int ncons = a_divF.nComp();
  a_divF.setVal(0.0);
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Real bndryArea = m_ebisBox.bndryArea(vof);
      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          Real bndryFlux = a_ebIrregFlux(vof,ivar);
          Real update = 0.;
          for ( int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& fluxDir = a_centroidFlux[idir];
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  Vector<FaceIndex> faces =
                    m_ebisBox.getFaces(vof, idir, sit());
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& face = faces[iface];
                      Real areaFrac = m_ebisBox.areaFrac(face);
                      Real faceFlux = fluxDir(face, ivar);
                      update += isign*areaFrac*faceFlux/m_dx[idir];
                    }
                  ch_flops()+=faces.size()*4; 
                }
            }

          //add EB boundary conditions in divergence
          update += bndryFlux*bndryArea*m_dx[0];
         //note NOT divided by volfrac
          a_divF(vof, ivar) = update;
        } //end loop over variables
    } //end loop over vofs
}


/******************/
///for when eb irreg flux is always zero
void
EBAdvectPatchIntegrator::
consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const IntVectSet&      a_ivs)
{
  CH_TIME("EBADvectPatchIntegrator::consUndividedDivergence");
  int ncons = a_divF.nComp();
  a_divF.setVal(0.0);
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          Real update = 0.;
          for ( int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& fluxDir = a_centroidFlux[idir];
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  Vector<FaceIndex> faces =
                    m_ebisBox.getFaces(vof, idir, sit());
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& face = faces[iface];
                      Real areaFrac = m_ebisBox.areaFrac(face);
                      Real faceFlux = fluxDir(face, ivar);
                      update += isign*areaFrac*faceFlux/m_dx[idir];
                    }
                  ch_flops()+=faces.size()*4; 
                }
            }

         //note NOT divided by volfrac
          a_divF(vof, ivar) = update;
        } //end loop over variables
    } //end loop over vofs
}
/******************/
FaceStencil
EBAdvectPatchIntegrator::
getInterpStencil(const FaceIndex& a_face) const
{
  
  FaceStencil sten;
  Real areaFrac = m_ebisBox.areaFrac(a_face);
  Real fullEnoughFaceAreaFrac = 0.9999999;
  //for full faces, no interpolation is necessary.
  if(areaFrac > fullEnoughFaceAreaFrac)
    {
      sten.add(a_face, 1.);
    }
  else
    {
      sten = EBArith::getInterpStencil(a_face, m_cfivs, 
                                       m_ebisBox, m_domain.domainBox());
    }
  return sten;
}
/************/
void
EBAdvectPatchIntegrator::
interpolateFluxToCentroids(BaseIFFAB<Real>              a_centroidFlux[SpaceDim],
                           const BaseIFFAB<Real>* const a_fluxInterpolant[SpaceDim],
                           const IntVectSet&            a_irregIVS)
{
  CH_TIME("EBAdvectPatchIntegrator::interpolateFluxToCentroids");
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

  //could just use numFluxes but this allows unit testing
  //with a single variable interpolant
  int nflux = a_centroidFlux[0].nComp();
  //now loop through the irregular faces
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      BaseIFFAB<Real>& fluxDir = a_centroidFlux[faceDir];
      const BaseIFFAB<Real>& interpol = *(a_fluxInterpolant[faceDir]);
      for (FaceIterator faceit(a_irregIVS, m_ebisBox.getEBGraph(), faceDir, stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& centFace = faceit();
          FaceStencil sten = getInterpStencil(centFace);
          for (int ivar = 0; ivar < nflux; ivar++)
            {
              Real newFlux = 0.0;
              for (int isten = 0; isten < sten.size(); isten++)
                {
                  const FaceIndex& stenFace= sten.face(isten);
                  Real weight = sten.weight(isten);
                  Real interpFlux = interpol(stenFace, ivar);
                  newFlux += weight*interpFlux;
                }
              fluxDir(centFace, ivar) = newFlux;
            }
          ch_flops()+=sten.size()*nflux*2; 
        }
    }
}

/*******************/
void
EBAdvectPatchIntegrator::
extrapolatePrim(EBFluxFAB&                       a_flux,
                Vector<BaseIVFAB<Real> * >&      a_coveredFluxMinu,
                Vector<BaseIVFAB<Real> * >&      a_coveredFluxPlus,
                const Vector<IntVectSet >&       a_coveredSetsMinu,
                const Vector<IntVectSet >&       a_coveredSetsPlus,
                const Vector<Vector<VolIndex> >& a_coveredFaceMinu,
                const Vector<Vector<VolIndex> >& a_coveredFacePlus,
                const EBCellFAB&                 a_consState,
                const EBCellFAB&                 a_source,
                const DataIndex&                 a_dit,
                const Real&                      a_time,
                const Real&                      a_dt)
{
  CH_TIME("EBAdvectPatchIntegrator::extrapolatePrim");

  EBCellFAB primMinu[SpaceDim];
  EBCellFAB primPlus[SpaceDim];
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      primMinu[idir].define(m_ebisBox, m_validBoxG4, 1);
      primPlus[idir].define(m_ebisBox, m_validBoxG4, 1);
    }

#if CH_SPACEDIM==2

  extrapolatePrim2D(primMinu, primPlus,
                    a_consState, a_source,
                    a_dit, a_time, a_dt);

#elif CH_SPACEDIM==3

  extrapolatePrim3D(primMinu   ,   primPlus,
                    a_consState,  a_source,
                    a_dit, a_time, a_dt);

#else
  MayDay::Error("bogus SpaceDim");
#endif

  //this keeps the fluxes from being calculated
  //on boundaries
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      Box bndryFaceBox;
      Box faceBox;
      bndryFaceBox = m_validBox;
      bndryFaceBox &= m_domain;
      bndryFaceBox.surroundingNodes(idir);

      faceBox = m_validBox;
      faceBox.grow(idir,1);
      faceBox &= m_domain;
      faceBox.grow(idir,-1);
      faceBox.surroundingNodes(idir);

      //define extended state
      BaseIVFAB<Real> extendStateMinuG4(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), 1);
      BaseIVFAB<Real> extendStatePlusG4(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), 1);
      extendStatePlusG4.setVal(0.);
      extendStateMinuG4.setVal(0.);

      //extrapolate to covered faces using updated extrapolated state
      extrapToCoveredFaces(extendStateMinuG4,
                           primMinu[idir],
                           primPlus[idir],
                           a_consState,
                           a_coveredFaceMinu[idir],
                           idir, Side::Lo, m_validBox);

      extrapToCoveredFaces(extendStatePlusG4,
                           primMinu[idir],
                           primPlus[idir],
                           a_consState,
                           a_coveredFacePlus[idir],
                           idir, Side::Hi, m_validBox);


      riemann(a_flux[idir], primPlus[idir], primMinu[idir], idir, faceBox);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(a_flux, a_consState, primMinu[idir],
                   Side::Lo,  a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox, idir);
      m_bc->fluxBC(a_flux, a_consState, primPlus[idir],
                   Side::Hi,  a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox, idir);
      //solve riemann problem between extended state and
      //value in vof to get covered flux.
      riemann(*a_coveredFluxMinu[idir],
              extendStateMinuG4, primMinu[idir],
              a_coveredFaceMinu[idir], idir, Side::Lo, m_validBox);
      riemann(*a_coveredFluxPlus[idir],
              extendStatePlusG4, primPlus[idir],
              a_coveredFacePlus[idir], idir, Side::Hi, m_validBox);
    }
}
/******/
void
EBAdvectPatchIntegrator::
extrapolatePrim2D(EBCellFAB           a_primMinu[SpaceDim],
                  EBCellFAB           a_primPlus[SpaceDim],
                  const EBCellFAB&    a_consState,
                  const EBCellFAB&    a_source,
                  const DataIndex&    a_dit,
                  const Real&         a_time,
                  const Real&         a_dt)
{
  CH_TIME("EBADvectPatchIntegrator::extrapolatePrim2D");

  //now define the plethora of data holders that I need
  
  EBFaceFAB        fluxOne[SpaceDim];
  BaseIVFAB<Real>  coveredFluxNormMinu[SpaceDim];
  BaseIVFAB<Real>  coveredFluxNormPlus[SpaceDim];
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      fluxOne[idir].define(m_ebisBox, m_validBoxG4, idir, 1);
      coveredFluxNormPlus[idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), 1);
      coveredFluxNormMinu[idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), 1);
    }

  //2.  compute slopes delta^d w
  //compute normal derivative stuff

  doNormalDerivativeExtr2D(a_primMinu,
                           a_primPlus,
                           fluxOne,
                           coveredFluxNormMinu,
                           coveredFluxNormPlus,
                           m_coveredFaceMinuG4,
                           m_coveredFacePlusG4,
                           a_consState,
                           a_source,
                           a_dit,  a_time, a_dt);


  /**/
  // Do the final corrections to the fluxes
  //this only happens on the non-ghosted box
  finalExtrap2D(a_primMinu,
                a_primPlus,
                coveredFluxNormMinu,
                coveredFluxNormPlus,
                m_coveredFaceMinuG4,
                m_coveredFacePlusG4,
                fluxOne,
                a_consState, a_time, a_dt);

  /**/



}
/*******************/
void
EBAdvectPatchIntegrator::
doNormalDerivativeExtr2D(EBCellFAB              a_primMinu[SpaceDim],
                         EBCellFAB              a_primPlus[SpaceDim],
                         EBFaceFAB              a_fluxOne[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                         const EBCellFAB&       a_primState,
                         const EBCellFAB&       a_source,
                         const DataIndex&       a_dit,
                         const Real     &       a_time,
                         const Real     &       a_dt)

{
  CH_TIME("EBAPI::doNormalDerivative2D");
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      Box slopeBoxG1 = grow(m_validBox, 1);
      Box slopeBoxG2 = grow(m_validBox, 2);
      slopeBoxG1 &= m_domain;
      slopeBoxG2 &= m_domain;
      Box cellBoxG2;
      Box cellBoxG1;
      Box faceBox;
      Box bndryFaceBox;
      int numSlop = 1;

      cellBoxG2 = m_validBox;
      cellBoxG2.grow(2);
      cellBoxG2.grow(idir, -1);
      cellBoxG2 &= m_domain;

      cellBoxG1 = m_validBox;
      cellBoxG1.grow(1);
      cellBoxG1.grow(idir, -1);
      cellBoxG1 &= m_domain;

      bndryFaceBox = m_validBox;
      bndryFaceBox.grow(1);
      bndryFaceBox &= m_domain;
      bndryFaceBox.surroundingNodes(idir);


      faceBox = m_validBox;
      faceBox.grow(2);
      faceBox &= m_domain;
      faceBox.grow(idir, -1);
      faceBox.surroundingNodes(idir);

      EBCellFAB slopesPrim(m_ebisBox, slopeBoxG2, numSlop);

      BaseIVFAB<Real> extendStateNormPlus(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), 1);
      BaseIVFAB<Real> extendStateNormMinu(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), 1);

      extendStateNormPlus.setVal(0.);
      extendStateNormMinu.setVal(0.);

      a_coveredFluxNormPlus[idir].setVal(0.);
      a_coveredFluxNormMinu[idir].setVal(0.);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
      a_fluxOne [idir].setVal(0.);

      slope(slopesPrim, a_primState, idir, slopeBoxG2);


      //3. Compute the effect of the normal derivative terms and the
      //source term on the extrapolation in space and time from cell
      //centers to faces. Equation 1.7 (top).
      normalPred(a_primMinu[idir], a_primPlus[idir],
                 a_primState, slopesPrim,
                 a_dt/m_dx[idir], idir, slopeBoxG2);

      // If the source term is valid add it to the primitive quantities
      //if source term is not defined, assume it is zero.  lots of
      //problems will have no sources.  not the most elegant method
      //but it beats passing a null pointer.
      if (a_source.isDefined())
        {
          incrementWithSource(a_primMinu[idir], a_source,
                              0.5*a_dt, slopeBoxG2);
          incrementWithSource(a_primPlus[idir], a_source,
                              0.5*a_dt, slopeBoxG2);
        }

      //4. compute estimates of the flux suitable for computing
      // 1D flux derivatives using a riemann solver for the interior R
      //and for the boundary RB.  Equation 1.8
      riemann(a_fluxOne[idir], a_primPlus[idir], a_primMinu[idir], idir, faceBox);

      //some wackiness to work around the fact that fluxfab did
      //not exist when a lot of this was written.
      Vector<EBFaceFAB*> bcfluxes(SpaceDim);
      for (int kdir = 0; kdir < SpaceDim; kdir++)
        {
          bcfluxes[kdir] = &(a_fluxOne[kdir]);
        }
      EBFluxFAB fluxAlias;
      fluxAlias.alias(bcfluxes);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(fluxAlias, a_primState, a_primMinu[idir],
                   Side::Lo, a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox, idir);
      m_bc->fluxBC(fluxAlias, a_primState, a_primPlus[idir],
                   Side::Hi, a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox, idir);

      //extrapolate to covered faces
      extrapToCoveredFaces(extendStateNormMinu,
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormMinu[idir],
                           idir, Side::Lo, slopeBoxG1);

      extrapToCoveredFaces(extendStateNormPlus,
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormPlus[idir],
                           idir, Side::Hi, slopeBoxG1);

      //solve riemann problem between extended state and
      //value in vof to get covered flux  Equation 1.9
      riemann(a_coveredFluxNormMinu[idir],
                extendStateNormMinu      , a_primMinu[idir],
              a_coveredFaceNormMinu[idir], idir, Side::Lo, slopeBoxG2);
      riemann(a_coveredFluxNormPlus[idir],
                extendStateNormPlus      , a_primPlus[idir],
              a_coveredFaceNormPlus[idir], idir, Side::Hi, slopeBoxG2);
    }
}
/*****************************/
void
EBAdvectPatchIntegrator::
slope(EBCellFAB&       a_slopePrim,
      const EBCellFAB& a_primState,
      const int&       a_dir,
      const Box&       a_box)
{
  CH_TIMERS("EBADvectPatchIntegrator::slope");
  CH_TIMER("advect_slopes:fortran",t1);
  CH_TIMER("advect_slopes:irregular",t2);

  //merging stuff from EBADvectPatchIntegrator::doSecondOrderSlopes
  int numSlope;
  Box box1, box2;
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  numSlope = 1;
  box1 = a_box;
  box1.grow(a_dir, 1);
  box1 &= m_domain;

  box2 = a_box;
  box2.grow(a_dir, 2);
  box2 &= m_domain;

  eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box1, m_domain, a_dir);


  BaseFab<Real>&         regDelta2W = a_slopePrim.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();

  int iuseLimiting = 0;
  if(m_useLimiting)
    {
      iuseLimiting = 1;
    }
  //this includes limiting if required (by m_useLimiting)
  CH_START(t1);
  FORT_ADVECTSLOPEDIFFS(CHF_FRA(regDelta2W),
                        CHF_CONST_FRA(regPrimState),
                        CHF_CONST_INT(numSlope),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(loBox),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(hiBox),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(centerBox),
                        CHF_INT(iuseLimiting));
  CH_STOP(t1);
  CH_START(t2);
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          int ivar = 0;

          bool hasFacesLeft, hasFacesRigh;
          Real dql, dqr, dqc;
          bool verbose =false;
          pointGetSlopes(dql, dqr,dqc,
                         hasFacesLeft,
                         hasFacesRigh,
                         vof, a_primState, a_dir, ivar, verbose);
          Real dqSec=0;
          if (!m_useLimiting)
            {
              dqSec = dqc;
            }
          else
            {
              if (hasFacesLeft && hasFacesRigh)
                {
                  Real dqlim = dqc;

                  FORT_POINTVLLIMITER(CHF_REAL(dqlim),
                                      CHF_CONST_REAL(dql),
                                      CHF_CONST_REAL(dqr));
                  dqSec = dqlim;
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  dqSec = 0.0;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  if (dqc*dql > 0.0)
                    {
                      Real rsign = 1.0;
                      if (dqc < 0.0)
                        {
                          rsign = -1.0;
                        }
                      dqSec = rsign*Min(Abs(dqc), Abs(dql));
                      ch_flops()+=5;
                    }
                  else
                    {
                      dqSec = 0.0;
                    }

                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  if (dqc*dqr > 0.0)
                    {
                      Real rsign = 1.0;
                      if (dqc < 0.0)
                        {
                          rsign = -1.0;
                        }
                      dqSec = rsign*Min(Abs(dqc), Abs(dqr));
                      ch_flops()+=4;
                    }
                  else
                    {
                      dqSec = 0.0;
                    }
                }
              else
                {
                  MayDay::Error("EBADvectPatchIntegrator::doSecondOrderSlopes -- missed a case");
                }
            }//end if using limiting

          a_slopePrim(vof, ivar) = dqSec;
        }
    }//end loop over vofs
  CH_STOP(t2);
}
/******************/
void
EBAdvectPatchIntegrator::
finalExtrap2D(EBCellFAB              a_primMinu[SpaceDim],
              EBCellFAB              a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
              const EBFaceFAB        a_fluxOne[SpaceDim],
              const EBCellFAB&       a_primState,
              const Real     &       a_time,
              const Real     &       a_dt)

{
  CH_TIME("EBADvectPatchIntegrator::finalExtrap2D");
  Box slopeBoxG1 = grow(m_validBox, 1);
  slopeBoxG1 &= m_domain;

  //7. compute fnal corrections to W due to the
  //final transverse derivatives.
  // Do the final corrections to the fluxes
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          // A different direction has been found
          if (diffDir != faceDir)
            {

              // In 2D,
              //the current primitive state is updated by a flux in
              // the other direction
              //a_coveredFlux holds the covered flux
              //due to the normal derivative
              //equation 1.18 line 3.
              updatePrim(a_primMinu[faceDir],
                         a_primPlus[faceDir], a_fluxOne[diffDir],
                         a_coveredFluxNormMinu[diffDir], a_coveredFluxNormPlus[diffDir],
                         a_coveredFaceNormMinu[diffDir], a_coveredFaceNormPlus[diffDir],
                         diffDir,  slopeBoxG1, (0.5)*a_dt/m_dx[diffDir]);

            } //end if dir2 != faceDir
        } //loop over dir2
    } //loop over facedir
}
///////////
void
EBAdvectPatchIntegrator::
extrapolatePrim3D(EBCellFAB           a_primMinu[SpaceDim],
                  EBCellFAB           a_primPlus[SpaceDim],
                  const EBCellFAB&    a_consState,
                  const EBCellFAB&    a_source,
                  const DataIndex&    a_dit,
                  const Real&         a_time,
                  const Real&         a_dt)
{
  CH_TIME("EBADvectPatchIntegrator::extrapolatePrim3D");
  //define the plethora of data holders that I need
  EBFaceFAB        fluxOne[SpaceDim];
  EBFaceFAB        fluxTwo[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  coveredFluxNormMinu[SpaceDim];
  BaseIVFAB<Real>  coveredFluxNormPlus[SpaceDim];
  BaseIVFAB<Real>  coveredFluxMinu3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  coveredFluxPlus3D[SpaceDim][SpaceDim];
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      fluxOne[idir].define(m_ebisBox, m_validBoxG4, idir, 1);
      coveredFluxNormPlus[idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), 1);
      coveredFluxNormMinu[idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), 1);
      for(int jdir = 0; jdir < SpaceDim; jdir++)
        {
          fluxTwo[idir][jdir].define(m_ebisBox, m_validBoxG4, idir, 1);
          if (idir != jdir) //the diagonal ones do not get used.
            {
              coveredFluxPlus3D[idir][jdir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), 1);
              coveredFluxMinu3D[idir][jdir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), 1);
            }
        }
    }


  //2.  compute slopes delta^d w
  //compute normal derivative stuff
  doNormalDerivativeExtr3D(a_primMinu,
                           a_primPlus,
                           fluxOne,
                           coveredFluxNormMinu,
                           coveredFluxNormPlus,
                           m_coveredFaceMinuG4,
                           m_coveredFacePlusG4,
                           a_consState,
                           a_source,
                           a_dit, a_time, a_dt);

  //6. In 3D compute corrections to u corresponding to one set
  //   of transverse derivatives appropriate to obtain (1, 1, 1)
  //   coupling.
  // In 3D, compute some additional intermediate fluxes
  // NOTE:  The diagonal entries of this array of fluxes are not
  // used and will not be defined.

  do111coupling(fluxTwo,
                coveredFluxMinu3D,
                coveredFluxPlus3D,
                a_primMinu,
                a_primPlus,
                coveredFluxNormMinu,
                coveredFluxNormPlus,
                m_coveredFaceMinuG4,
                m_coveredFacePlusG4,
                fluxOne,
                a_consState,
                a_dit, a_time, a_dt);

  finalExtrap3D(a_primMinu,
                a_primPlus,
                coveredFluxMinu3D,
                coveredFluxPlus3D,
                fluxTwo,
                a_consState, a_time, a_dt);

}
/*******************/
void 
EBAdvectPatchIntegrator::
doNormalDerivativeExtr3D(EBCellFAB              a_primMinu[SpaceDim],
                         EBCellFAB              a_primPlus[SpaceDim],
                         EBFaceFAB              a_fluxOne[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                         const EBCellFAB&       a_primState,
                         const EBCellFAB&       a_source,
                         const DataIndex&       a_dit,
                         const Real     &       a_time,
                         const Real     &       a_dt)
{
  CH_TIME("EBADvectPatchIntegrator::doNormalDerivative3D");
  int numSlop = 1;

  //all one big loop since no cross terms here
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      Box modBoxCov;
      Box modBoxOpen;
      Box faceBox;
      Box bndryFaceBox;

      //set up the data structures
      modBoxOpen = m_validBox;
      modBoxOpen.grow(3);
      modBoxOpen.grow(idir, -1);
      modBoxOpen &= m_domain;

      modBoxCov = m_validBox;
      modBoxCov.grow(2);
      modBoxCov.grow(idir, -1);
      modBoxCov &= m_domain;

      bndryFaceBox = modBoxCov;
      bndryFaceBox.surroundingNodes(idir);

      faceBox = modBoxCov;
      faceBox.grow(idir, 1);
      faceBox &= m_domain;
      faceBox.grow(idir, -1);
      faceBox.surroundingNodes(idir);

      BaseIVFAB<Real> extendStateNormPlus(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), 1);
      BaseIVFAB<Real> extendStateNormMinu(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), 1);

      extendStateNormPlus.setVal(0.);
      extendStateNormMinu.setVal(0.);

      extendStateNormPlus.setVal(0.);
      extendStateNormMinu.setVal(0.);
      a_coveredFluxNormPlus[idir].setVal(0.);
      a_coveredFluxNormMinu[idir].setVal(0.);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
      a_fluxOne[idir].setVal(0.);

      //compute slopes
      EBCellFAB slopesPrim(m_ebisBox, modBoxOpen, numSlop);
      slope(slopesPrim, a_primState, idir, modBoxOpen);

      //3. Compute the effect of the normal derivative terms and the
      //source term on the extrapolation in space and time from cell
      //centers to faces. Equation 1.7 (top).
      normalPred(a_primMinu[idir], a_primPlus[idir],
                 a_primState, slopesPrim,
                 a_dt/m_dx[idir], idir, modBoxOpen);

      // If the source term is valid add it to the primitive quantities
      //if source term is not defined, assume it is zero.  lots of
      //problems will have no sources.  not the most elegant method
      //but it beats passing a null pointer.
      if (a_source.isDefined())
        {
          incrementWithSource(a_primMinu[idir], a_source,
                              0.5*a_dt, modBoxOpen);
          incrementWithSource(a_primPlus[idir], a_source,
                              0.5*a_dt, modBoxOpen);
        }

      //4. compute estimates of the flux suitable for computing
      // 1D flux derivatives using a riemann solver for the interior R
      //and for the boundary RB.  Equation 1.8
      riemann(a_fluxOne[idir], a_primPlus[idir], a_primMinu[idir],
              idir, faceBox);

      //some hackiness to work around the fact that fluxfab did
      //not exist when a lot of this was written.
      Vector<EBFaceFAB*> bcfluxes(SpaceDim);
      for (int kdir = 0; kdir < SpaceDim; kdir++)
        {
          bcfluxes[kdir] = &(a_fluxOne[kdir]);
        }
      EBFluxFAB fluxAlias;
      fluxAlias.alias(bcfluxes);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(fluxAlias, a_primState, a_primMinu[idir],
                   Side::Lo, a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox, idir);
      m_bc->fluxBC(fluxAlias, a_primState, a_primPlus[idir],
                   Side::Hi, a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox, idir);

      //extrapolate to covered faces
      extrapToCoveredFaces(extendStateNormMinu,
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormMinu[idir],
                           idir, Side::Lo, modBoxCov);

      extrapToCoveredFaces(extendStateNormPlus,
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormPlus[idir],
                           idir, Side::Hi, modBoxCov);

      //solve riemann problem between extended state and
      //value in vof to get covered flux  Equation 1.9
      riemann(a_coveredFluxNormMinu[idir],
              extendStateNormMinu, a_primMinu[idir],
              a_coveredFaceNormMinu[idir], idir, Side::Lo, modBoxOpen);
      riemann(a_coveredFluxNormPlus[idir],
              extendStateNormPlus, a_primPlus[idir],
              a_coveredFaceNormPlus[idir], idir, Side::Hi, modBoxOpen);
    }
}

/*******************/
void 
EBAdvectPatchIntegrator::
do111coupling(EBFaceFAB              a_fluxTwo[SpaceDim][SpaceDim],
              BaseIVFAB<Real>        a_coveredFluxMinu3D[SpaceDim][SpaceDim],
              BaseIVFAB<Real>        a_coveredFluxPlus3D[SpaceDim][SpaceDim],
              const EBCellFAB        a_primMinu[SpaceDim],
              const EBCellFAB        a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
              const EBFaceFAB        a_fluxOne[SpaceDim],
              const EBCellFAB&       a_primState,
              const DataIndex&       a_dit,
              const Real     &       a_time,
              const Real     &       a_dt)
{
  CH_TIME("EBADvectPatchIntegrator::do111coupling");
  Box slopeBoxG1 = grow(m_validBox, 1);
  Box slopeBoxG2 = grow(m_validBox, 2);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
  Box faceBox[SpaceDim][SpaceDim];
  Box bndryFaceBox[SpaceDim][SpaceDim];
  Box modBoxCov[SpaceDim];
  Box modBoxOpen[SpaceDim][SpaceDim];
  for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
    {
      for (int diffDir = 0; diffDir < SpaceDim; ++diffDir)
        {

          modBoxOpen[faceDir][diffDir] = grow(m_validBox, 2);
          modBoxOpen[faceDir][diffDir].grow(diffDir, -1);
          modBoxOpen[faceDir][diffDir] &= m_domain;
        }
    }

  for (int modDir = 0; modDir < SpaceDim; ++modDir)
    {
      modBoxCov[modDir] = grow(m_validBox, 1);
      modBoxCov[modDir] &= m_domain;
      for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
        {
          bndryFaceBox[faceDir][modDir] = modBoxCov[modDir];
          bndryFaceBox[faceDir][modDir].surroundingNodes(faceDir);

          faceBox[faceDir][modDir] = modBoxCov[modDir];
          faceBox[faceDir][modDir].grow(faceDir, 1);
          faceBox[faceDir][modDir] &= m_domain;
          faceBox[faceDir][modDir].grow(faceDir, -1);
          faceBox[faceDir][modDir].surroundingNodes(faceDir);
        }
    }

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          if (diffDir != faceDir)
            {
              a_coveredFluxPlus3D[faceDir][diffDir].setVal(0.);
              a_coveredFluxMinu3D[faceDir][diffDir].setVal(0.);
              a_fluxTwo[faceDir][diffDir].setVal(0.);

              BaseIVFAB<Real> extendStatePlus3D(m_coveredSetsPlusG4[faceDir], m_ebisBox.getEBGraph(), 1);
              BaseIVFAB<Real> extendStateMinu3D(m_coveredSetsMinuG4[faceDir], m_ebisBox.getEBGraph(), 1);
              extendStateMinu3D.setVal(0.);
              extendStatePlus3D.setVal(0.);

              //we are modify the third direction
              int modDir = 3 - faceDir - diffDir;

              // Copy data for in place modification
              EBCellFAB primMinuTemp, primPlusTemp;
              primMinuTemp.clone(a_primMinu[faceDir]);
              primPlusTemp.clone(a_primPlus[faceDir]);

              // Update the current, extrapolated primitive
              //variable using a flux
              // in a different direction.  uses covered fluxes calculated
              // above.  Equation 1.12
              updatePrim(primMinuTemp, primPlusTemp, a_fluxOne[diffDir],
                         a_coveredFluxNormMinu[diffDir], a_coveredFluxNormPlus[diffDir],
                         a_coveredFaceNormMinu[diffDir], a_coveredFaceNormPlus[diffDir],
                         diffDir,  modBoxOpen[faceDir][diffDir],  (1.0/3.0)*a_dt/m_dx[diffDir]);

              extrapToCoveredFaces(extendStateMinu3D,
                                   primMinuTemp,
                                   primPlusTemp,
                                   a_primState,
                                   m_coveredFaceMinuG4[faceDir],
                                   faceDir, Side::Lo, modBoxCov[modDir]);

              extrapToCoveredFaces(extendStatePlus3D,
                                   primMinuTemp,
                                   primPlusTemp,
                                   a_primState,
                                   m_coveredFacePlusG4[faceDir],
                                   faceDir, Side::Hi, modBoxCov[modDir]);

              // Solve the Riemann problem and get fluxes.  Eqution 1.15
              riemann(a_fluxTwo[faceDir][diffDir], primPlusTemp, primMinuTemp,
                      faceDir, faceBox[faceDir][modDir]);

              //some wackiness to work around the fact that fluxfab did
              //not exist when a lot of this was written.
              Vector<EBFaceFAB*> bcfluxes(SpaceDim);
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  bcfluxes[idir] = &(a_fluxTwo[idir][diffDir]);
                }
              EBFluxFAB fluxAlias;
              fluxAlias.alias(bcfluxes);

              m_bc->fluxBC(fluxAlias,  a_primState, primMinuTemp,
                           Side::Lo,  a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox[faceDir][diffDir], faceDir);
              m_bc->fluxBC(fluxAlias,  a_primState, primPlusTemp,
                           Side::Hi,  a_time, m_ebisBox, a_dit,m_validBox, bndryFaceBox[faceDir][diffDir], faceDir);

              //solve riemann problem between extended state and
              //value in vof to get covered flux.  Equation 1.14
              riemann(a_coveredFluxMinu3D[faceDir][diffDir],
                      extendStateMinu3D, primMinuTemp,
                      m_coveredFaceMinuG4[faceDir], faceDir, Side::Lo, modBoxOpen[faceDir][diffDir]);
              riemann(a_coveredFluxPlus3D[faceDir][diffDir],
                      extendStatePlus3D, primPlusTemp,
                      m_coveredFacePlusG4[faceDir], faceDir, Side::Hi, modBoxOpen[faceDir][diffDir]);

            }
        }
    }
}
/***********/
void 
EBAdvectPatchIntegrator::
finalExtrap3D(EBCellFAB              a_primMinu[SpaceDim],
              EBCellFAB              a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxMinu3D[SpaceDim][SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxPlus3D[SpaceDim][SpaceDim],
              const EBFaceFAB        a_fluxTwo[SpaceDim][SpaceDim],
              const EBCellFAB&       a_primState,
              const Real     &       a_time,
              const Real     &       a_dt)

{
  CH_TIME("EBADvectPatchIntegrator::finalExtrap3D");
  Box slopeBoxG1 = grow(m_validBox, 1);
  slopeBoxG1 &= m_domain;

  //7. compute fnal corrections to W due to the
  //final transverse derivatives.
  // Do the final corrections to the fluxes
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
          // A different direction has been found
          if (dir2 != dir1)
            {
              // In 3D,  find a direction different from the two above
              // here we need to use the coveredFlux1D for the covered flux

              int dir3 = 3 - dir2  - dir1;
              // A different direction has been found
              //const EBFluxFAB& fluxTwo = a_fluxTwo[dir3];
              // Update the conservative state
              //using both corrected fluxes in
              // the other two directions.  Equation 1.18 line 4
              updatePrim(a_primMinu[dir1], a_primPlus[dir1],
                         a_fluxTwo[dir2][dir3],
                         a_coveredFluxMinu3D[dir2][dir3],
                         a_coveredFluxPlus3D[dir2][dir3],
                         m_coveredFaceMinuG4[dir2],
                         m_coveredFacePlusG4[dir2],
                         dir2,  slopeBoxG1,  (0.5)*a_dt/m_dx[dir2]);

            }
        } //loop over dir2

    } //loop over dir1
}
///begin working set
/*******************/
void
EBAdvectPatchIntegrator::
extrapolateBCG(EBFluxFAB&                       a_flux,
               const EBCellFAB&                 a_consState,
               const EBCellFAB&                 a_source,
               const DataIndex&                 a_dit,
               const Real&                      a_time,
               const Real&                      a_dt)
{
  CH_TIMERS("EBAdvectPatchIntegrator::extrapolateBCG(fluxfab)");
  CH_TIMER("define_left_right_primitive_vars",t1);
  CH_TIMER("taylor_extrapolation",t2);
  CH_TIMER("solve_riemann",t3);

  CH_START(t1);
  EBCellFAB primMinu[SpaceDim];
  EBCellFAB primPlus[SpaceDim];
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      primMinu[idir].define(m_ebisBox, m_validBoxG4, 1);
      primPlus[idir].define(m_ebisBox, m_validBoxG4, 1);
    }
  CH_STOP(t1);

  CH_START(t2);
  //Taylor extrapolation
  extrapolateBCG(primMinu,
                 primPlus,
                 a_consState,
                 a_source,
                 a_dit, a_time, a_dt);
  CH_STOP(t2);

  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  //this keeps the fluxes from being calculated on boundaries
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      bndryFaceBox[idir] = m_validBox;
      bndryFaceBox[idir] &= m_domain;
      bndryFaceBox[idir].surroundingNodes(idir);
      faceBox[idir] = m_validBox;
      faceBox[idir].grow(idir,1);
      faceBox[idir] &= m_domain;
      faceBox[idir].grow(idir,-1);
      faceBox[idir].surroundingNodes(idir);
    }

  CH_START(t3);
  //solve the Riemann problem for plus and minus states at a face
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      riemann(a_flux[idir], primPlus[idir], primMinu[idir],
              idir, faceBox[idir]);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      // m_bc->fluxBC(a_flux, a_primState, m_primMinu[idir],
      m_bc->fluxBC(a_flux, a_consState, primMinu[idir],
                   Side::Lo,  a_time, m_ebisBox, a_dit, m_validBox,
                   bndryFaceBox[idir], idir);
      // m_bc->fluxBC(a_flux, a_primState, m_primPlus[idir],
      m_bc->fluxBC(a_flux, a_consState, primPlus[idir],
                   Side::Hi,  a_time, m_ebisBox, a_dit, m_validBox,
                   bndryFaceBox[idir], idir);
    }
  CH_STOP(t3);

}

/******/
void
EBAdvectPatchIntegrator::
extrapolateBCG(EBCellFAB           a_primMinu[SpaceDim],
               EBCellFAB           a_primPlus[SpaceDim],
               const EBCellFAB&    a_consState,
               const EBCellFAB&    a_source,
               const DataIndex&    a_dit,
               const Real&                      a_time,
               const Real&                      a_dt)
{
  CH_TIMERS("EBADvectPatchIntegrator::extrapolateBCG(cellfab)");
  CH_TIMER("slopesNorm",t1);
  CH_TIMER("incrementWithSource",t2);
  CH_TIMER("slopesTang",t3);
  CH_TIMER("floorPrimitives",t4);

  Box slopeBoxG1 = grow(m_validBox, 1);
  Box slopeBoxG2 = grow(m_validBox, 2);
  Box slopeBoxG3 = grow(m_validBox, 3);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
  slopeBoxG3 &= m_domain;
  Box cellBoxG2[SpaceDim];
  Box cellBoxG1[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  int numSlop = 1;
  EBCellFAB slopesPrim[SpaceDim];
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      cellBoxG2[idir] = m_validBox;
      cellBoxG2[idir].grow(2);
      cellBoxG2[idir].grow(idir, -1);
      cellBoxG2[idir] &= m_domain;

      cellBoxG1[idir] = m_validBox;
      cellBoxG1[idir].grow(1);
      cellBoxG1[idir].grow(idir, -1);
      cellBoxG1[idir] &= m_domain;

      bndryFaceBox[idir] = m_validBox;
      bndryFaceBox[idir].grow(1);
      bndryFaceBox[idir] &= m_domain;
      bndryFaceBox[idir].surroundingNodes(idir);

      faceBox[idir] = m_validBox;
      faceBox[idir].grow(2);
      faceBox[idir] &= m_domain;
      faceBox[idir].grow(idir, -1);
      faceBox[idir].surroundingNodes(idir);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      slopesPrim[idir].define(m_ebisBox, slopeBoxG1, numSlop);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      EBCellFAB normSlopes;
      EBCellFAB tanSlopes[SpaceDim-1];

      Box normSlopeBox;
      Box tanSlopeBox[SpaceDim-1];
      normSlopeBox = m_validBox;
      normSlopeBox.grow(idir, 1);
      normSlopeBox &= m_domain;

      CH_START(t1);
      normSlopes.define(m_ebisBox, normSlopeBox, numSlop);
      slope(normSlopes, a_consState, idir, normSlopeBox);
      normalPred(a_primMinu[idir],
                 a_primPlus[idir],
                 a_consState,
                 normSlopes,
                 a_dt/m_dx[idir],
                 idir,
                 normSlopeBox);
      CH_STOP(t1);

      CH_START(t2);
      if (a_source.isDefined())
        {
          incrementWithSource(a_primMinu[idir],
                              a_source,
                              0.5*a_dt,
                              normSlopeBox);
          incrementWithSource(a_primPlus[idir],
                              a_source,
                              0.5*a_dt,
                              normSlopeBox);
        }
      CH_STOP(t2);

      CH_START(t3);
      Tuple<int,SpaceDim-1> tanDirs = PolyGeom::computeTanDirs(idir);
      for (int itan = 0; itan < SpaceDim-1; itan++)
        {
          int tanDir = tanDirs[itan];
          tanSlopeBox[itan] = normSlopeBox;
          tanSlopeBox[itan].grow(tanDir, 1);
          tanSlopeBox[itan] &= m_domain;

          tanSlopes[itan].define(m_ebisBox, tanSlopeBox[itan], numSlop);

          //Minion stability fix
          EBCellFAB statePlusSource;
          statePlusSource.clone(a_consState);
          statePlusSource.plus(a_source, 0.5*a_dt);

          upwindSlope(tanSlopes[itan],
                      statePlusSource,
                      tanDir,
                      tanSlopeBox[itan]);

          transversePred(a_primMinu[idir],
                         a_primPlus[idir],
                         a_consState,
                         tanSlopes[itan],
                         a_dt/m_dx[tanDir],
                         tanDir,
                         normSlopeBox);
        }
      CH_STOP(t3);

      CH_START(t4);
      floorPrimitives(a_primMinu[idir], normSlopeBox);
      floorPrimitives(a_primPlus[idir], normSlopeBox);
      CH_STOP(t4);

    }
}
///end working set
/**************
the box here is grown from m_valid box
***************/
void
EBAdvectPatchIntegrator::
averageVelToCC(EBCellFAB&                        a_normalVel,
               const EBFluxFAB&                  a_advectionVel,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloLo,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloHi,
               const Vector<Vector<VolIndex> >&  a_coveredFaceLo,
               const Vector<Vector<VolIndex> >&  a_coveredFaceHi,
               const Box&                        a_box) const
{
  CH_TIME("EBAdvectPatchIntegrator::averageVelToCC");
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const EBFaceFAB& faceNormVel = a_advectionVel[idir];

      //treat every cell as regular
      BaseFab<Real>&       regCellVel = a_normalVel.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceNormVel.getSingleValuedFAB();
      FORT_EBAVEFACETOCELL(CHF_FRA(regCellVel),
                         CHF_CONST_FRA1(regFaceVel,0),
                         CHF_CONST_INT(idir),
                         CHF_BOX(a_box));
    }
  //correct way
  {
    CH_TIME("EBADvectPatchIntegrator::irregularLoop");
    for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        //the set can be bigger than the box for performance reasons.
        if (a_box.contains(vof.gridIndex()))
          {
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, idir, Side::Lo);
                Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, idir, Side::Hi);
                Real velLo= 0.0;
                if (facesLo.size() > 0)
                  {
                    for (int iface = 0; iface < facesLo.size(); iface++)
                      {
                        velLo += a_advectionVel[idir](facesLo[iface], 0);
                      }
                    velLo /=  facesLo.size();
                  }
                else
                  {
                    const BaseIVFAB<Real>& coveredVelLo = *a_coveredVeloLo[idir];
                    velLo = coveredVelLo(vof, 0);
                  }
                Real velHi= 0.0;
                if (facesHi.size() > 0)
                  {
                    for (int iface = 0; iface < facesHi.size(); iface++)
                      {
                        velHi += a_advectionVel[idir](facesHi[iface], 0);
                      }
                    velHi /=  facesHi.size();
                  }
                else
                  {
                    const BaseIVFAB<Real>& coveredVelHi = *a_coveredVeloHi[idir];
                    velHi = coveredVelHi(vof, 0);
                  }

                Real velAve = 0.5*(velHi + velLo);
                a_normalVel(vof, idir) = velAve;
              }//end loop over directions
          }
      }//end loop over irregular vofs
  }
}
/******/
void
EBAdvectPatchIntegrator::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBFaceFAB&        a_primFace,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBAdvectPatchIntegrator::extrapToCoveredFaces");
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = 1;
          Vector<Real> extPrim(numPrim, 0.0);
#if CH_SPACEDIM==2
          pointExtrapToCovered2D(extPrim,
                                 a_primFace,
                                 a_primState,
                                 a_faceDir,
                                 vof,
                                 normal,
                                 a_sd,
                                 numPrim);
#elif CH_SPACEDIM==3
          pointExtrapToCovered3D(extPrim,
                                 a_primFace,
                                 a_primState,
                                 a_faceDir,
                                 vof,
                                 normal,
                                 a_sd,
                                 numPrim);
#else
          bogus_ch_spacedim();
#endif
          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              a_extendedPrim(vof, ivar) = extPrim[ivar];
            }
        }
    }
  const IntVectSet& ivs  = a_extendedPrim.getIVS();
  floorPrimitives(a_extendedPrim, ivs);
}

/********/
void
EBAdvectPatchIntegrator::
pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                       const EBFaceFAB&        a_primFace,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  CH_assert(SpaceDim==2);
  int tangenDir = 1 - a_faceDir;

  int signNorm = 1;
  int signTang = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  if (Abs(a_normal[tangenDir]) == 0.0)
    {//since signTang is arbitrary for this case,
      //check if we are next to the high domain edge
      IntVect ivHi = a_vof.gridIndex();;
      ivHi[tangenDir] += 1;
      if (!m_domain.contains(ivHi))
        {
          signTang = -1;
        }
    }
  else if (a_normal[tangenDir] < 0.0)
    {
      signTang = -1;
    }

  //iv[0][0] is the vof. iv[1][0] is vofside.  iv[0][1] is vofup iv[1][1] is vofcorner
  IntVect   ivSten[2][2];
  bool      hasVoF[2][2];
  VolIndex     vof[2][2];
  FaceIndex   face[2][2];
  Real         val[2][2];
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      for (int ix = 0; ix < 2; ix++)
        {
          for (int iy = 0; iy < 2; iy++)
            {
              Vector<FaceIndex> facesTmp;
              hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], vofsStencil, ivSten[ix][iy]);
              if (hasVoF[ix][iy])
                {
                  facesTmp = m_ebisBox.getFaces(vof[ix][iy], a_faceDir, flip(a_sd));
                  hasVoF[ix][iy] = (hasVoF[ix][iy] && (facesTmp.size() == 1));
                }
              if (hasVoF[ix][iy])
                {
                  face[ix][iy] = facesTmp[0];
                  val [ix][iy] = a_primFace(face[ix][iy], ivar);
                }
              else
                {
                  val[ix][iy] = 0.0;
                }
            }
        }

      int dirBigNorm, dirLitNorm;
      VolIndex vofBigNorm, vofLitNorm;
      bool hasVoFBigNorm,  hasVoFLitNorm;
      int  signBigNorm, signLitNorm;
      Real wBigNorm, wLitNorm;
      Real dWBigNorm[2],dWLitNorm[2];
      if (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir])
        {
          dirBigNorm = tangenDir;
          dirLitNorm = a_faceDir;

          signBigNorm   = signTang;
          signLitNorm   = signNorm;
          hasVoFBigNorm = hasVoF[0][1];
          if (hasVoF[0][1])
            {
              wBigNorm      =    val[0][1];
              vofBigNorm    =    vof[0][1];
            }
          else
            {
              wBigNorm      =    val[0][0];
              vofBigNorm    =    vof[0][0];
            }
        }
      else
        {
          dirBigNorm = a_faceDir;
          dirLitNorm = tangenDir;

          signBigNorm   = signNorm;
          signLitNorm   = signTang;
          hasVoFBigNorm = hasVoF[1][0];

          if (hasVoF[1][0])
            {
              //              Real deltaW = a_slopesSeco[a_faceDir](vof[1][0]  , ivar);
              Real deltaW;
              coveredExtrapSlopes(deltaW, vof[1][0],a_primState, a_faceDir, ivar);
              wBigNorm      =    val[1][0] - signBigNorm*deltaW;
              vofBigNorm    =    vof[1][0];
            }
          else
            {
              wBigNorm    =    val[0][0];
              vofBigNorm  =    vof[0][0];
            }

        }
      hasVoFLitNorm = hasVoF[1][1];
      if (hasVoF[1][1])
        {
          wLitNorm      =    val[1][1];
          vofLitNorm    =    vof[1][1];
        }
      else
        {
          wLitNorm      =    val[0][0];
          vofLitNorm    =    vof[0][0];
        }

      if (hasVoFLitNorm && hasVoFBigNorm)
        {
          coveredExtrapSlopes(dWBigNorm[0],vofBigNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWBigNorm[1],vofBigNorm,a_primState,1,ivar);
          coveredExtrapSlopes(dWLitNorm[0],vofLitNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWLitNorm[1],vofLitNorm,a_primState,1,ivar);

          const Real xLit = Abs( (m_dx[dirBigNorm]*a_normal[dirLitNorm])/ (m_dx[dirLitNorm]*a_normal[dirBigNorm]) );
          const Real xBig = 1.0 - xLit;

          const Real dWBigExtrap = xLit*dWLitNorm[dirBigNorm] + xBig*dWBigNorm[dirBigNorm];
          const Real dWLitExtrap = xLit*dWLitNorm[dirLitNorm] + xBig*dWBigNorm[dirLitNorm];

          const Real wc = xLit*wLitNorm + xBig*wBigNorm;
          a_extrapVal[ivar] = wc -signBigNorm*dWBigExtrap - xLit*signLitNorm*dWLitExtrap;
          // a_extrapVal[ivar] = wc;
        }
      else if (hasVoF[1][0])
        {//change methods, still second order
          Real deltaW;
          coveredExtrapSlopes(deltaW,vof[1][0],a_primState,a_faceDir,ivar);

          a_extrapVal[ivar] = val[1][0] - 2.0*signNorm*deltaW;
          // a_extrapVal[ivar] = val[1][0];
        }
      else if (hasVoFLitNorm)
        {//drop order
          //pout() << "EBAdvectPatchIntegrator::pointExtrapToCovered2D:lit: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = wLitNorm;
        }
      else if (hasVoFBigNorm)
        {
          //pout() << "EBAdvectPatchIntegrator::pointExtrapToCovered2D:big: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = wBigNorm;
        }
      else
        {
          //pout() << "EBAdvectPatchIntegrator::pointExtrapToCovered2D:none: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = a_primState(a_vof,ivar);
        }
    }
}

void
EBAdvectPatchIntegrator::
pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                       const EBFaceFAB&        a_primFace,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  RealVect normal = a_normal;
  //HACK put in true normal
  //RealVect truenorm = getTrueNorm(a_vof,a_faceDir,a_sd);
  //normal = truenorm;
  //END HACK

  CH_assert(SpaceDim==3);
  Tuple<int,CH_SPACEDIM-1> tangenDir = PolyGeom::computeTanDirs(a_faceDir);
  Real anormNorm = Abs(normal[a_faceDir]);
  Real anormTan[2];
  int signNorm = 1;
  if (normal[a_faceDir] < 0.0) signNorm = -1;
  int signTang[2];
  for (int itan = 0; itan < 2; itan++)
    {
      int tanDir = tangenDir[itan];
      anormTan[itan] = Abs(normal[tanDir]);
      signTang[itan] =  1;
      if (anormTan[itan] == 0.0)
        {//since signTang is arbitrary for this case,
          //check if we are next to the high domain edge
          IntVect ivHi = a_vof.gridIndex();;
          ivHi[tanDir] += 1;
          if (!m_domain.contains(ivHi))
            {
              signTang[itan] = -1;
            }
        }
      else if (normal[tanDir] < 0.0)
        {
          signTang[itan] = -1;
        }
    }

  const IntVect& ivVoF= a_vof.gridIndex();

  //whice one of the tangential directions has the largest normal
  int d1, d2;
  if (anormTan[0]/m_dx[tangenDir[0]] > anormTan[1]/m_dx[tangenDir[1]])
    {
      d1 = 0;
      d2 = 1;
    }
  else
    {
      d1 = 1;
      d2 = 0;
    }

  // figure out in which plane we are extrapolating
  bool faceDirOut = ((anormNorm/m_dx[a_faceDir] > anormTan[0]/m_dx[tangenDir[0]]) &&
                     (anormNorm/m_dx[a_faceDir] > anormTan[1]/m_dx[tangenDir[1]]));

  IntVect ivSten[2][2];
  bool hasVoF[2][2];
  VolIndex    vofSten[2][2];
  FaceIndex  faceSten[2][2];

  if (faceDirOut)
    {// face direction has largest normal.

      ivSten[0][0] = ivVoF + signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[0]*BASISV(tangenDir[0]);
      ivSten[0][1] = ivVoF + signTang[1]*BASISV(tangenDir[1]);
      ivSten[1][1] = ivVoF + signTang[0]*BASISV(tangenDir[0]) + signTang[1]*BASISV(tangenDir[1]) ;
      d1 = 0;
      d2 = 1;

    }
  else
    { //tandir[d1] is the biggest

      ivSten[0][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      - signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      ;
      ivSten[0][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) - signNorm*BASISV(a_faceDir);
      ivSten[1][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) ;

    }

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  //get the face info for a 1D extrap stencil (if needed)
  VolIndex vofSten1D;
  FaceIndex faceSten1D;
  const IntVect ivSten1D = ivVoF + signNorm*BASISV(a_faceDir);
  bool has1DSten = EBArith::isVoFHere(vofSten1D, vofsStencil, ivSten1D);
  Vector<FaceIndex> facesTmp1D;
  if (has1DSten)
    {
      facesTmp1D = m_ebisBox.getFaces(vofSten1D, a_faceDir, flip(a_sd));
      if (facesTmp1D.size() == 1)
        {
          faceSten1D = facesTmp1D[0];
        }
      else
        {
          has1DSten = false;
        }
    }

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], vofsStencil, ivSten[ix][iy]);
          Vector<FaceIndex> facesTmp;
          if (hasVoF[ix][iy])
            {
             facesTmp = m_ebisBox.getFaces(vofSten[ix][iy], a_faceDir, flip(a_sd));
             hasVoF[ix][iy] = (hasVoF[ix][iy] && (facesTmp.size() == 1));
            }
          if (hasVoF[ix][iy])
            {
              faceSten[ix][iy] = facesTmp[0];
            }
        }
    }
  bool hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      if (hasAllVoFs)
        {
          Real WVal[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  WVal[ix][iy] = a_primFace(faceSten[ix][iy], ivar);
                }
            }
          Real delWNorm[2][2];
          Real delWTan0[2][2];
          Real delWTan1[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2 ; iy++)
                {
                  coveredExtrapSlopes(delWNorm[ix][iy], vofSten[ix][iy], a_primState, a_faceDir   , ivar);
                  coveredExtrapSlopes(delWTan0[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d1], ivar);
                  coveredExtrapSlopes(delWTan1[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d2], ivar);
                }
            }
          if (faceDirOut)
            {
              // case 00 is special because we have to move it backwards in face direction to get it to plane
              Real deltaW;
              coveredExtrapSlopes(deltaW, vofSten[0][0],a_primState, a_faceDir, ivar);
              //Real slipSlope = a_slopesSeco[a_faceDir](vofSten[0][0], ivar);
              Real slipSlope = deltaW;
              WVal[0][0] -= signNorm*slipSlope;
              //debug to make compatible with 2d
              //WVal[0][1] = WVal[0][0];
              //end debug

              Real x0Val = (m_dx[a_faceDir]*anormTan[0])/(m_dx[tangenDir[0]]*anormNorm);
              Real x1Val = (m_dx[a_faceDir]*anormTan[1])/(m_dx[tangenDir[1]]*anormNorm);
              Real    funcVal = bilinearFunc(    WVal,x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);
              Real dWdn =  -signNorm*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*x0Val*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

            }
          else
            {
              //tanDir[d1] the biggest normal
              Real x0Val = (m_dx[tangenDir[d1]]*anormNorm   )/(m_dx[    a_faceDir]*anormTan[d1]);
              Real x1Val = (m_dx[tangenDir[d1]]*anormTan[d2])/(m_dx[tangenDir[d2]]*anormTan[d1]);
              Real    funcVal = bilinearFunc(WVal,    x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);

              Real dWdn =  -signNorm*x0Val*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

            }
        }
      else if (has1DSten)
        {//try 1D extrap
          const Real WVal1D = a_primFace(faceSten1D, ivar);
          Real deltaW;
          coveredExtrapSlopes(deltaW,vofSten1D,a_primState,a_faceDir,ivar);
          a_extrapVal[ivar] = WVal1D - 2.0*signNorm*deltaW;
        }
      else
        // At least one of the vofs in the stencil exists - we drop order
        // by using a weighted sum of the the wVal's.
        {
          //pout() << "EBADvectPatchIntegrator::pointExtrapToCovered3D: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = 0.;
          Real volTot = 0.;
          bool hasAVoF = false;
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  if (hasVoF[ix][iy])
                    {
                      hasAVoF = true;
                      Real wVal = a_primFace(faceSten[ix][iy], ivar);
                      Real volFrac = m_ebisBox.volFrac(vofSten[ix][iy]);
                      a_extrapVal[ivar] += wVal*volFrac;
                      volTot += volFrac;
                    }
                }
            }
          if (hasAVoF)
            {
              CH_assert(volTot > 0.0);
              a_extrapVal[ivar] = a_extrapVal[ivar]/volTot;
            }
          else
            {
              a_extrapVal[ivar] = a_primState(a_vof, ivar);
            }
        }
    }
}
/******/
void
EBAdvectPatchIntegrator::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBCellFAB&        a_primMinu,
                     const EBCellFAB&        a_primPlus,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBADvectPatchIntegrator::extrapToCoveredFaces");
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = 1;
          Vector<Real> extPrim(numPrim, 0.0);
          if (SpaceDim== 2)
            {
              pointExtrapToCovered2D(extPrim,
                                     a_primMinu,
                                     a_primPlus,
                                     a_primState,
                                     a_faceDir,
                                     vof,
                                     normal,
                                     a_sd,
                                     numPrim);
            }
          else if (SpaceDim==3)
            {
              pointExtrapToCovered3D(extPrim,
                                     a_primMinu,
                                     a_primPlus,
                                     a_primState,
                                     a_faceDir,
                                     vof,
                                     normal,
                                     a_sd,
                                     numPrim);
            }
          else
            {
              MayDay::Error("Bogus SpaceDim");
            }

          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              a_extendedPrim(vof, ivar) = extPrim[ivar];
            }
        }
    }

  const IntVectSet& ivs  = a_extendedPrim.getIVS();
  floorPrimitives(a_extendedPrim, ivs);
}
/********/
void
EBAdvectPatchIntegrator::
pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  CH_assert(SpaceDim==2);
  int tangenDir = 1 - a_faceDir;

  int signNorm = 1;
  int signTang = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  if (Abs(a_normal[tangenDir]) == 0.0)
    {//since signTang is arbitrary for this case,
      //check if we are next to the high domain edge
      IntVect ivHi = a_vof.gridIndex();;
      ivHi[tangenDir] += 1;
      if (!m_domain.contains(ivHi))
        {
          signTang = -1;
        }
    }
  else if (a_normal[tangenDir] < 0.0)
    {
      signTang = -1;
    }

  //iv[0][0] is the vof. iv[1][0] is vofside.  iv[0][1] is vofup iv[1][1] is vofcorner
  IntVect   ivSten[2][2];
  bool      hasVoF[2][2];
  VolIndex     vof[2][2];
  Real         val[2][2];
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      for (int ix = 0; ix < 2; ix++)
        {
          for (int iy = 0; iy < 2; iy++)
            {
              hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], vofsStencil, ivSten[ix][iy]);
              if (hasVoF[ix][iy])
                {
                  if (a_sd == Side::Hi)
                    {
                      val[ix][iy] = a_primMinu(vof[ix][iy], ivar);
                    }
                  else
                    {
                      val[ix][iy] = a_primPlus(vof[ix][iy], ivar);
                    }
                }
              else
                {
                  val[ix][iy] = 0.0;
                }
            }
        }

      int dirBigNorm, dirLitNorm;
      VolIndex vofBigNorm, vofLitNorm;
      bool hasVoFBigNorm,  hasVoFLitNorm;
      int  signBigNorm, signLitNorm;
      Real wBigNorm, wLitNorm;
      Real dWBigNorm[2],dWLitNorm[2];
      Real eps = 1.0e-12;
      if ((Abs(a_normal[a_faceDir]) < eps) || (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir]))
        {
          dirBigNorm = tangenDir;
          dirLitNorm = a_faceDir;

          signBigNorm   = signTang;
          signLitNorm   = signNorm;
          hasVoFBigNorm = hasVoF[0][1];
          if (hasVoF[0][1])
            {
              wBigNorm      =    val[0][1];
              vofBigNorm    =    vof[0][1];
            }
          else
            {
              wBigNorm      =    val[0][0];
              vofBigNorm    =    vof[0][0];
            }
        }
      else
        {
          dirBigNorm = a_faceDir;
          dirLitNorm = tangenDir;

          signBigNorm   = signNorm;
          signLitNorm   = signTang;
          hasVoFBigNorm = hasVoF[1][0];

          if (hasVoF[1][0])
            {
              //              Real deltaW = a_slopesSeco[a_faceDir](vof[1][0]  , ivar);
              Real deltaW;
              coveredExtrapSlopes(deltaW, vof[1][0],a_primState, a_faceDir, ivar);
              wBigNorm      =    val[1][0] - signBigNorm*deltaW;
              vofBigNorm    =    vof[1][0];
            }
          else
            {
              wBigNorm    =    val[0][0];
              vofBigNorm  =    vof[0][0];
            }

        }
      hasVoFLitNorm = hasVoF[1][1];
      if (hasVoF[1][1])
        {
          wLitNorm      =    val[1][1];
          vofLitNorm    =    vof[1][1];
        }
      else
        {
          wLitNorm      =    val[0][0];
          vofLitNorm    =    vof[0][0];
        }

      if (hasVoFLitNorm && hasVoFBigNorm)
        {
          coveredExtrapSlopes(dWBigNorm[0],vofBigNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWBigNorm[1],vofBigNorm,a_primState,1,ivar);
          coveredExtrapSlopes(dWLitNorm[0],vofLitNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWLitNorm[1],vofLitNorm,a_primState,1,ivar);

          const Real xLit = Abs( (m_dx[dirBigNorm]*a_normal[dirLitNorm])/ (m_dx[dirLitNorm]*a_normal[dirBigNorm]) );
          const Real xBig = 1.0 - xLit;

          const Real dWBigExtrap = xLit*dWLitNorm[dirBigNorm] + xBig*dWBigNorm[dirBigNorm];
          const Real dWLitExtrap = xLit*dWLitNorm[dirLitNorm] + xBig*dWBigNorm[dirLitNorm];

          const Real wc = xLit*wLitNorm + xBig*wBigNorm;
          a_extrapVal[ivar] = wc - xBig*signBigNorm*dWBigExtrap - xLit*signLitNorm*dWLitExtrap;
          //debug turn off extrapolation bit
          //a_extrapVal[ivar] = wc;
          //end debug
          //no flattening here
        }
      else if (hasVoF[1][0])
        {//change methods, still second order
          Real deltaW;
          coveredExtrapSlopes(deltaW,vof[1][0],a_primState,a_faceDir,ivar);

          a_extrapVal[ivar] = val[1][0] - 2.0*signNorm*deltaW;
        }
      else if (hasVoFLitNorm)
        {//drop order
          a_extrapVal[ivar] = wLitNorm;
        }
      else if (hasVoFBigNorm)
        {
          a_extrapVal[ivar] = wBigNorm;
        }
      else
        {
          a_extrapVal[ivar] = a_primState(a_vof,ivar);
        }
    }
}
/*******/
void
EBAdvectPatchIntegrator::
pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  RealVect normal = a_normal;
  CH_assert(SpaceDim==3);
  Tuple<int,CH_SPACEDIM-1> tangenDir = PolyGeom::computeTanDirs(a_faceDir);
  Real anormNorm = Abs(normal[a_faceDir]);
  Real anormTan[2];
  int signNorm = 1;
  if (normal[a_faceDir] < 0.0) signNorm = -1;
  int signTang[2];
  for (int itan = 0; itan < 2; itan++)
    {
      int tanDir = tangenDir[itan];
      anormTan[itan] = Abs(normal[tanDir]);
      signTang[itan] =  1;
      if (anormTan[itan] == 0.0)
        {//since signTang is arbitrary for this case,
          //check if we are next to the high domain edge
          IntVect ivHi = a_vof.gridIndex();;
          ivHi[tanDir] += 1;
          if (!m_domain.contains(ivHi))
            {
              signTang[itan] = -1;
            }
        }
      else if (normal[tanDir] < 0.0)
        {
          signTang[itan] = -1;
        }
    }

  const IntVect& ivVoF= a_vof.gridIndex();

  //whice one of the tangential directions has the largest normal
  int d1, d2;
  if (anormTan[0]/m_dx[tangenDir[0]] > anormTan[1]/m_dx[tangenDir[1]])
    {
      d1 = 0;
      d2 = 1;
    }
  else
    {
      d1 = 1;
      d2 = 0;
    }

  // figure out in which plane we are extrapolating
  bool faceDirOut = ((anormNorm/m_dx[a_faceDir] > anormTan[0]/m_dx[tangenDir[0]]) &&
                     (anormNorm/m_dx[a_faceDir] > anormTan[1]/m_dx[tangenDir[1]]));

  IntVect ivSten[2][2];
  bool hasVoF[2][2];
  VolIndex vofSten[2][2];

  if (faceDirOut)
    {// face direction has largest normal.

      ivSten[0][0] = ivVoF + signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[0]*BASISV(tangenDir[0]);
      ivSten[0][1] = ivVoF + signTang[1]*BASISV(tangenDir[1]);
      ivSten[1][1] = ivVoF + signTang[0]*BASISV(tangenDir[0]) + signTang[1]*BASISV(tangenDir[1]) ;
      d1 = 0;
      d2 = 1;

    }
  else
    { //tandir[d1] is the biggest

      ivSten[0][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      - signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      ;
      ivSten[0][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) - signNorm*BASISV(a_faceDir);
      ivSten[1][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) ;

    }

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], vofsStencil, ivSten[ix][iy]);
        }
    }
  bool hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  //get the vof info for a 1D extrap stencil (if needed)
  VolIndex vofSten1D;
  const IntVect ivSten1D = ivVoF + signNorm*BASISV(a_faceDir);
  const bool has1DVoF =  EBArith::isVoFHere(vofSten1D, vofsStencil, ivSten1D);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      if (hasAllVoFs)
        {
          Real WVal[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  if (a_sd == Side::Hi)
                    {
                      WVal[ix][iy] = a_primMinu(vofSten[ix][iy], ivar);
                    }
                  else
                    {
                      WVal[ix][iy] = a_primPlus(vofSten[ix][iy], ivar);
                    }
                }

            }

          Real delWNorm[2][2];
          Real delWTan0[2][2];
          Real delWTan1[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2 ; iy++)
                {
                  coveredExtrapSlopes(delWNorm[ix][iy], vofSten[ix][iy], a_primState, a_faceDir   , ivar);
                  coveredExtrapSlopes(delWTan0[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d1], ivar);
                  coveredExtrapSlopes(delWTan1[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d2], ivar);
                }
            }

          if (faceDirOut)
            {
              // case 00 is special because we have to move it backwards in face direction to get it to plane
              Real deltaW;
              coveredExtrapSlopes(deltaW, vofSten[0][0],a_primState, a_faceDir, ivar);
              //Real slipSlope = a_slopesSeco[a_faceDir](vofSten[0][0], ivar);
              Real slipSlope = deltaW;
              WVal[0][0] -= signNorm*slipSlope;
              //debug to make compatible with 2d
              //WVal[0][1] = WVal[0][0];
              //end debug

              Real x0Val = (m_dx[a_faceDir]*anormTan[0])/(m_dx[tangenDir[0]]*anormNorm);
              Real x1Val = (m_dx[a_faceDir]*anormTan[1])/(m_dx[tangenDir[1]]*anormNorm);
              Real    funcVal = bilinearFunc(    WVal,x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);
              Real dWdn =  -signNorm*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*x0Val*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;
              //no flattening here
            }
          else
            {
              //tanDir[d1] the biggestnormal
              Real x0Val = (m_dx[tangenDir[d1]]*anormNorm   )/(m_dx[    a_faceDir]*anormTan[d1]);
              Real x1Val = (m_dx[tangenDir[d1]]*anormTan[d2])/(m_dx[tangenDir[d2]]*anormTan[d1]);
              Real    funcVal = bilinearFunc(WVal,    x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);

              Real dWdn =  -signNorm*x0Val*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

              //no flattening here
            }
        }
      else if (has1DVoF)
        {//try 1D extrap
          Real WVal1D;
          if (a_sd == Side::Hi)
            {
              WVal1D = a_primMinu(vofSten1D, ivar);
            }
          else
            {
              WVal1D = a_primPlus(vofSten1D, ivar);
            }
          Real deltaW;
          coveredExtrapSlopes(deltaW,vofSten1D,a_primState,a_faceDir,ivar);
          a_extrapVal[ivar] = WVal1D - 2.0*signNorm*deltaW;
        }
      else
        // At least one of the vofs in the stencil exists - we drop order
        // by using a weighted sum of the the wVal's.
        {
          a_extrapVal[ivar] = 0.;
          Real volTot = 0.;
          bool hasAVoF = false;
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  if (hasVoF[ix][iy])
                    {
                      hasAVoF = true;
                      Real wVal;
                      if (a_sd == Side::Hi)
                        {
                          wVal = a_primMinu(vofSten[ix][iy], ivar);
                        }
                      else
                        {
                          wVal = a_primPlus(vofSten[ix][iy], ivar);
                        }
                      Real volFrac = m_ebisBox.volFrac(vofSten[ix][iy]);
                      a_extrapVal[ivar] += wVal*volFrac;
                      volTot += volFrac;
                    }
                }
            }
          if (hasAVoF)
            {
              CH_assert(volTot > 0.0);
              a_extrapVal[ivar] = a_extrapVal[ivar]/volTot;
            }
          else
            {
              a_extrapVal[ivar] = a_primState(a_vof, ivar);
            }
        }
      // None of the vofs in the stencil exists. We use the value at the
      // cell center adjacent to the face.
    }
  //    if (!hasAllVoFs)
  //          cout << "Dropping order at covered face = " << a_vof <<
  //          "direction = " << a_faceDir <<  "\n";
}
/*****************************/
void 
EBAdvectPatchIntegrator::
coveredExtrapSlopes(Real&            a_dq,
                    const VolIndex&  a_vof,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const int&       a_ivar)
{
  CH_TIME("EBADvectPatchIntegrator::coveredExtrapSlopes");
  bool hasFacesLeft, hasFacesRigh;
  Real dql, dqr;
  bool verbose = false;
  pointGetSlopes(dql, dqr,a_dq,
                 hasFacesLeft,
                 hasFacesRigh,
                 a_vof, a_primState, a_dir, a_ivar, verbose);
  //no flattening here.
}
/*************Need this because don't want to limit transverse slopes in BCG****************/
void 
EBAdvectPatchIntegrator::
pointGetSlopesUpwind(Real&            a_dql,
                     Real&            a_dqr,
                     Real&            a_dqc,
                     bool&            a_hasFacesLeft,
                     bool&            a_hasFacesRigh,
                     const VolIndex&  a_vof,
                     const EBCellFAB& a_primState,
                     const int&       a_dir,
                     const int&       a_ivar,
                     const bool&      a_verbose)
{
  CH_TIME("EBPG::pointGetSlopes");
  const IntVect&  iv = a_vof.gridIndex();
  //one-sided diffs on domain bndry
  const Box& domainBox = m_domain.domainBox();
  // const EBCellFAB& velcc = *m_normalVelPtr;
  bool onLeftDomain = (iv[a_dir] == domainBox.smallEnd(a_dir));
  bool onRighDomain = (iv[a_dir] == domainBox.bigEnd(a_dir)  );
  VolIndex vofLeft;
  VolIndex vofRigh;
  a_hasFacesLeft = (m_ebisBox.numFaces(a_vof, a_dir, Side::Lo)==1) && !onLeftDomain;
  a_hasFacesRigh = (m_ebisBox.numFaces(a_vof, a_dir, Side::Hi)==1) && !onRighDomain;

  if (a_hasFacesLeft)
    {
      Vector<FaceIndex> facesLeft =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Lo);
      vofLeft = facesLeft[0].getVoF(Side::Lo);
    }

  if (a_hasFacesRigh)
    {
      Vector<FaceIndex> facesRigh =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Hi);
      vofRigh = facesRigh[0].getVoF(Side::Hi);
    }

  Real valCent = a_primState(a_vof, a_ivar);
  Real valLeft = 0.0;
  Real valRigh = 0.0;
  if (a_hasFacesLeft)
    {
      valLeft = a_primState(vofLeft, a_ivar);
    }
  if (a_hasFacesRigh)
    {
      valRigh = a_primState(vofRigh, a_ivar);
    }

  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.
  //if neither exists, set all slopes to zero
  if (a_hasFacesLeft)
    {
      a_dql = valCent - valLeft;
      ch_flops()+=1; 
    }
  if (a_hasFacesRigh)
    {
      a_dqr = valRigh - valCent;
      ch_flops()+=1; 
    }
  if (a_hasFacesLeft && a_hasFacesRigh)
    {
      a_dqc = 0.5*(a_dql+a_dqr);
      ch_flops()+=2; 
    }
  else if (!a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dql = 0.0;
      a_dqr = 0.0;
      a_dqc = 0.0;
    }
  else if (a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dqr = a_dql;
      //no higher-order one sided diffs
      a_dqc = a_dql;
    }
  else if (a_hasFacesRigh && !a_hasFacesLeft)
    {
      a_dql = a_dqr;
      //no higher-order one sided diffs
      a_dqc = a_dqr;
    }
  else
    {
      MayDay::Error("EBADvectPatchIntegrator::pointGetSlopes -- missed a case");
    }
  if (a_verbose)
    {
      pout() << "  a_vof="    << a_vof
             << ", vofLeft= " << vofLeft
             << ", vofRigh= " << vofRigh  << endl;
      pout() << "  hasFacesLeft=" << a_hasFacesLeft
             << ", hasFacesRigh=" << a_hasFacesRigh;
      pout() << ", valLeft=" << valLeft
             << ", valRigh=" << valRigh
             << ", valCent=" << valCent << endl;
      pout() << "  a_dql=" << a_dql
             << ", a_dqr=" << a_dqr
             << ", a_dqc=" << a_dqc << endl;
    }
  if (a_dqc != a_dqc)
    {
      MayDay::Error("EBADvectPatchIntegrator::pointGetSlopesUpwind -- a_dqc != a_dqc");
    }
}
/*****************************/
void 
EBAdvectPatchIntegrator::
pointGetSlopes(Real&            a_dql,
               Real&            a_dqr,
               Real&            a_dqc,
               bool&            a_hasFacesLeft,
               bool&            a_hasFacesRigh,
               const VolIndex&  a_vof,
               const EBCellFAB& a_primState,
               const int&       a_dir,
               const int&       a_ivar,
               const bool&      a_verbose)
{
   const IntVect&  iv = a_vof.gridIndex();
  //one-sided diffs on domain bndry
  const Box& domainBox = m_domain.domainBox();
  bool onLeftDomain = (iv[a_dir] == domainBox.smallEnd(a_dir));
  bool onRighDomain = (iv[a_dir] == domainBox.bigEnd(a_dir)  );
  VolIndex vofLeft;
  VolIndex vofRigh;
  a_hasFacesLeft = (m_ebisBox.numFaces(a_vof, a_dir, Side::Lo)==1) && !onLeftDomain;
  a_hasFacesRigh = (m_ebisBox.numFaces(a_vof, a_dir, Side::Hi)==1) && !onRighDomain;

  if (a_hasFacesLeft)
    {
      Vector<FaceIndex> facesLeft =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Lo);
      vofLeft = facesLeft[0].getVoF(Side::Lo);
    }

  if (a_hasFacesRigh)
    {
      Vector<FaceIndex> facesRigh =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Hi);
      vofRigh = facesRigh[0].getVoF(Side::Hi);
    }

  Real valCent = a_primState(a_vof, a_ivar);
  Real valLeft = 0.0;
  Real valRigh = 0.0;
  if (a_hasFacesLeft)
    {
      valLeft = a_primState(vofLeft, a_ivar);
    }
  if (a_hasFacesRigh)
    {
      valRigh = a_primState(vofRigh, a_ivar);
    }

  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.
  //if neither exists, set all slopes to zero
  //limiting sprinked in here to make higher order
  //one-sided diffs possible
  if (a_hasFacesLeft)
    {
      a_dql = valCent - valLeft;
      ch_flops()+=1; 
    }
  if (a_hasFacesRigh)
    {
      a_dqr = valRigh - valCent;
      ch_flops()+=1; 
    }
  if (a_hasFacesLeft && a_hasFacesRigh)
    {
      a_dqc = 0.5*(a_dql+a_dqr);
      ch_flops()+=2; 
      if (m_useLimiting)
        {
          Real dqmin = Min(Abs(a_dql)*2.0,Abs(a_dqr)*2.0);
          if (a_dql*a_dqr < 0.)
            {
              dqmin = 0.;
              ch_flops()+=6; 
            }
          if (a_dqc < 0.)
            {
              a_dqc = -Min(dqmin,Abs(a_dqc));
              ch_flops()+=7; 
            }
          else
            {
              a_dqc = Min(dqmin,Abs(a_dqc));
              ch_flops()+=7; 
            }
        }
    }
  else if (!a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dql = 0.0;
      a_dqr = 0.0;
      a_dqc = 0.0;
    }
  else if (a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dqr = a_dql;
      //no higher-order one sided diffs
      a_dqc = a_dql;
    }
  else if (a_hasFacesRigh && !a_hasFacesLeft)
    {
      a_dql = a_dqr;
      //no higher-order one sided diffs
      a_dqc = a_dqr;
    }
  else
    {
      MayDay::Error("EBADvectPatchIntegrator::pointGetSlopes -- missed a case");
    }
  if (a_verbose)
    {
      pout() << "  a_vof="    << a_vof
             << ", vofLeft= " << vofLeft
             << ", vofRigh= " << vofRigh  << endl;
      pout() << "  hasFacesLeft=" << a_hasFacesLeft
             << ", hasFacesRigh=" << a_hasFacesRigh;
      pout() << ", valLeft=" << valLeft
             << ", valRigh=" << valRigh
             << ", valCent=" << valCent << endl;
      pout() << "  a_dql=" << a_dql
             << ", a_dqr=" << a_dqr
             << ", a_dqc=" << a_dqc << endl;
    }
  if (a_dqc != a_dqc)
    {
      MayDay::Error("EBADvectPatchIntegrator::pointGetSlopes -- a_dqc != a_dqc");
    }
  //no flattening here
}
  /// floors if m_isMaxMinSet
void
EBAdvectPatchIntegrator::
floorPrimitives(EBCellFAB&  a_primState,
                const Box&  a_box)
{
  if (m_isMaxMinSet)
    {
      // IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(a_box);
      IntVectSet ivsIrreg(a_box);
      for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          a_primState(vof, 0)  = Min(a_primState(vof, 0), m_maxVal);
          a_primState(vof, 0)  = Max(a_primState(vof, 0), m_minVal);
        }
    }
}

/// floors if m_isMaxMinSet
void
EBAdvectPatchIntegrator::
floorPrimitives(BaseIVFAB<Real>&   a_primState,
                const IntVectSet&  a_set)
{
  if (m_isMaxMinSet)
    {
      for (VoFIterator vofit(a_set, m_ebisBox.getEBGraph());vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          a_primState(vof, 0)  = Min(a_primState(vof, 0), m_maxVal);
          a_primState(vof, 0)  = Max(a_primState(vof, 0), m_minVal);
        }
    }
}
/******/
Real
EBAdvectPatchIntegrator::
bilinearFunc(const Real  a_WVal[2][2],
             const Real& a_xd1,
             const Real& a_xd2)
{
  Real D =                a_WVal[0][0];
  Real A =(a_WVal[1][0] - D);
  Real B =(a_WVal[0][1] - D);
  Real C =(a_WVal[1][1] - D - A - B);
  Real retval = A*a_xd1 + B*a_xd2 + C*a_xd1*a_xd2 + D;
  return retval;
}
///
/**
 */
void
EBAdvectPatchIntegrator::
normalPred(EBCellFAB&       a_rhoLo,
           EBCellFAB&       a_rhoHi,
           const EBCellFAB& a_rho,
           const EBCellFAB& a_dRho,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_TIME("EBADvectPatchIntegrator::normalPred");
  CH_assert(m_isVelSet);
  const EBCellFAB& velcc = *m_normalVelPtr;


  int ivar = 0;

  /**/
  BaseFab<Real>& regRhoLo    = a_rhoLo.getSingleValuedFAB();
  BaseFab<Real>& regRhoHi    = a_rhoHi.getSingleValuedFAB();
  const BaseFab<Real>& regVel =  velcc.getSingleValuedFAB();
  const BaseFab<Real>& regRho =  a_rho.getSingleValuedFAB();
  const BaseFab<Real>& regDRho= a_dRho.getSingleValuedFAB();

  FORT_PREDADVECT(CHF_BOX(a_box),
                  CHF_CONST_FRA1(regRho,ivar),
                  CHF_CONST_FRA1(regDRho, ivar),
                  CHF_CONST_FRA(regVel),
                  CHF_FRA1(regRhoLo, ivar),
                  CHF_FRA1(regRhoHi, ivar),
                  CHF_CONST_INT(a_dir),
                  CHF_CONST_REAL(a_dtbydx));
  /**/

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real dense, denlo, denhi, denslope;
      RealVect veloc;

      dense    = a_rho(vof, ivar);
      denslope = a_dRho(vof, ivar);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = velcc(vof, idir);
        }

      FORT_POINTPREDADVECT(CHF_REAL(dense),
                           CHF_REAL(denlo),
                           CHF_REAL(denhi),
                           CHF_REAL(denslope),
                           CHF_REALVECT(veloc),
                           CHF_CONST_INT(a_dir),
                           CHF_CONST_REAL(a_dtbydx));

      a_rhoLo(vof, ivar) = denlo;
      a_rhoHi(vof, ivar) = denhi;
    }
}

/******/
void
EBAdvectPatchIntegrator::
incrementWithSource(EBCellFAB&       a_state,
                    const EBCellFAB& a_source,
                    const Real&      a_scale,
                    const Box&       a_box)
{
  //if conservative source is set to true,  state is conservative state
  //if conservative source is set to false, state is primitive    state
  CH_TIME("EBADvectPatchIntegrator::incrmentWithSource");

  BaseFab<Real>&       primReg = a_state.getSingleValuedFAB();
  const BaseFab<Real>& sourReg =    a_source.getSingleValuedFAB();

  FORT_INCSOURCE(CHF_FRA(primReg),
                 CHF_CONST_FRA(sourReg),
                 CHF_CONST_REAL(a_scale),
                 CHF_BOX(a_box));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < a_source.nComp(); ivar++)
        {
          Real increm = a_scale*a_source(vof, ivar);
          a_state(vof, ivar) +=  increm;
        }
    }
}

/******/
void
EBAdvectPatchIntegrator::
updatePrim(EBCellFAB&              a_primMinu,
           EBCellFAB&              a_primPlus,
           const EBFaceFAB&        a_primFace,
           const BaseIVFAB<Real>&  a_coveredPrimMinu,
           const BaseIVFAB<Real>&  a_coveredPrimPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_faceDir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBADvectPatchIntegrator::updatePrim");
  //uses flux as prims
  //update the regular vofs.
  //first we have to copy the original state
  //so that the irregular stuff can be updated later.

  int numPrim = 1;
  //save state so that the regular can overwrite irregular cells
  EBCellFAB  cacheMinu, cachePlus;
  cacheMinu.clone(a_primMinu);
  cachePlus.clone(a_primPlus);

  const EBCellFAB& normalVel = *m_normalVelPtr;
  BaseFab<Real>& regPrimMinu = a_primMinu.getSingleValuedFAB();
  BaseFab<Real>& regPrimPlus = a_primPlus.getSingleValuedFAB();
  const BaseFab<Real>& regPrimFace      =    a_primFace.getSingleValuedFAB();

  {
    CH_TIME("EBADvectPatchIntegrator::fortran");
    const BaseFab<Real>& regNormVel = normalVel.getSingleValuedFAB();
    FORT_ADVECTUPDATE(CHF_BOX(a_box),
                      CHF_FRA(regPrimMinu),
                      CHF_FRA(regPrimPlus),
                      CHF_CONST_FRA(regPrimFace),
                      CHF_CONST_FRA(regNormVel),
                      CHF_CONST_INT(a_faceDir),
                      CHF_CONST_INT(numPrim),
                      CHF_CONST_REAL(a_scale),
                      CHF_BOX(a_box));
  }




  {
    CH_TIME("EBADvectPatchIntegrator::irregularLoop");
    for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        //the set can be bigger than the box for performance reasons.
        if (a_box.contains(vof.gridIndex()))
          {
            Real vel = normalVel(vof, a_faceDir);

            for (int ivar = 0; ivar < numPrim; ivar++)
              {

                Real primLo = 0;
                {
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, a_faceDir, Side::Lo);
                  if (facesLo.size() > 0)
                    {
                      for (int iface = 0; iface < facesLo.size(); iface++)
                        {
                          primLo += a_primFace(facesLo[iface], ivar);
                        }
                      primLo /= facesLo.size();
                    }
                  else
                    {
                      primLo = a_coveredPrimMinu(vof, ivar);
                    }
                }
                Real primHi = 0;
                {
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, a_faceDir, Side::Hi);
                  if (facesHi.size() > 0)
                    {
                      for (int iface = 0; iface < facesHi.size(); iface++)
                        {
                          primHi += a_primFace(facesHi[iface], ivar);
                        }
                      primHi /= facesHi.size();
                    }
                  else
                    {
                      primHi = a_coveredPrimPlus(vof, ivar);
                    }
                }
                Real primDiff = primHi - primLo;

                //state changed in regular update.
                //so we have to use the cache.
                //scale holds 0.5*dt/dx in 2d.
                a_primMinu(vof, ivar) =  cacheMinu(vof, ivar) - a_scale*vel*primDiff;
                a_primPlus(vof, ivar) =  cachePlus(vof, ivar) - a_scale*vel*primDiff;

              }
          }
      }
  }

  {
    CH_TIME("EBADvectPatchIntegrator::floors");
    floorPrimitives(a_primMinu, a_box);
    floorPrimitives(a_primPlus, a_box);
  }
}
void
EBAdvectPatchIntegrator::
upwindSlope(EBCellFAB&       a_slopeUpWi,
            const EBCellFAB& a_primState,
            const int&       a_dir,
            const Box&       a_box)
{
  CH_TIME("EBAPI::upwindSlope");

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               a_box, m_domain, a_dir);

  // a_slopeUpWi.define(m_ebisBox, entireBox, numSlope);

  BaseFab<Real>& regSlopeUpWi       = a_slopeUpWi.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
  const BaseFab<Real>& regNormalVel = m_normalVelPtr->getSingleValuedFAB();
  int numSlope = 1;

  FORT_UPWINDDIFFS(CHF_FRA(regSlopeUpWi),
                   CHF_CONST_FRA(regPrimState),
                   CHF_CONST_FRA(regNormalVel),
                   CHF_CONST_INT(numSlope),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          for (int ivar = 0; ivar < numSlope; ivar++)
            {
              bool hasFacesLeft, hasFacesRigh;
              Real dql, dqr, dqc;
              bool verbose =false;
              pointGetSlopesUpwind(dql, dqr, dqc,
                                   hasFacesLeft,
                                   hasFacesRigh,
                                   vof, a_primState, a_dir, ivar, verbose);


              Real velPt = (*m_normalVelPtr)(vof, a_dir);
              if (velPt > EBAPI_TOL)
                {
                  dqc = dql;
                }
              else if (velPt < -EBAPI_TOL)
                {
                  dqc = dqr;
                }
              else
                {
                  dqc = 0.5*(dql+dqr);
                  ch_flops()+=2;
                }

              a_slopeUpWi(vof, ivar) = dqc;
            } //end loop over variables
        }
    }//end loop over vofs
  {
    CH_TIME("bndry_slopes");
    //want to be able to call this in test codes
    if (m_isBCSet)
      {
        m_bc->setBndrySlopes(a_slopeUpWi, a_primState, m_ebisBox, entireBox, a_dir);
      }
  }
}
/*****************************/
void 
EBAdvectPatchIntegrator::
transversePred(EBCellFAB&       a_rhoLo,
               EBCellFAB&       a_rhoHi,
               const EBCellFAB& a_rho,
               const EBCellFAB& a_dRho,
               const Real&      a_dtbydx,
               const int&       a_dir,
               const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::transversePred");
  CH_assert(m_isVelSet);
  const EBCellFAB& velcc = *m_normalVelPtr;

  int ivar = 0;

  /**/
  BaseFab<Real>& regRhoLo    = a_rhoLo.getSingleValuedFAB();
  BaseFab<Real>& regRhoHi    = a_rhoHi.getSingleValuedFAB();
  const BaseFab<Real>& regVel =  velcc.getSingleValuedFAB();
  const BaseFab<Real>& regRho =  a_rho.getSingleValuedFAB();
  const BaseFab<Real>& regDRho= a_dRho.getSingleValuedFAB();

  FORT_PREDADVECTTRANS(CHF_BOX(a_box),
                       CHF_CONST_FRA1(regRho,ivar),
                       CHF_CONST_FRA1(regDRho, ivar),
                       CHF_CONST_FRA(regVel),
                       CHF_FRA1(regRhoLo, ivar),
                       CHF_FRA1(regRhoHi, ivar),
                       CHF_CONST_INT(a_dir),
                       CHF_CONST_REAL(a_dtbydx));
  /**/

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real dense, denlo, denhi, denslope;
      RealVect veloc;

      dense    = a_rho(vof, ivar);
      denslope = a_dRho(vof, ivar);
      denlo    = a_rhoLo(vof, ivar);
      denhi    = a_rhoHi(vof, ivar);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = velcc(vof, idir);
        }

      FORT_POINTPREDADVECTTRANS(CHF_REAL(dense),
                                CHF_REAL(denlo),
                                CHF_REAL(denhi),
                                CHF_REAL(denslope),
                                CHF_REALVECT(veloc),
                                CHF_CONST_INT(a_dir),
                                CHF_CONST_REAL(a_dtbydx));

      a_rhoLo(vof, ivar) = denlo;
      a_rhoHi(vof, ivar) = denhi;
    }
}

/****************/
void
EBAdvectPatchIntegrator::
advectiveDerivative(EBCellFAB&                       a_uDotDelRho,
                    const EBFluxFAB&                 a_faceRho,
                    const EBFluxFAB&                 a_faceVel,
                    const Vector<BaseIVFAB<Real>*> & a_coveredRhoLo,
                    const Vector<BaseIVFAB<Real>*> & a_coveredRhoHi,
                    const Vector<BaseIVFAB<Real>*> & a_coveredVelLo,
                    const Vector<BaseIVFAB<Real>*> & a_coveredVelHi,
                    const Vector<Vector<VolIndex> >& a_coveredFaceLo,
                    const Vector<Vector<VolIndex> >& a_coveredFaceHi,
                    const Box&                       a_box)
{
  CH_TIME("EBPatchAdvect::advectiveDerivative");
  int ncomp = a_faceRho.nComp();

  CH_assert(a_coveredRhoLo.size() == SpaceDim);
  CH_assert(a_coveredVelLo.size() == SpaceDim);
  CH_assert(a_coveredRhoHi.size() == SpaceDim);
  CH_assert(a_coveredVelHi.size() == SpaceDim);
  CH_assert(a_coveredFaceHi.size() == SpaceDim);
  CH_assert(a_coveredFaceLo.size() == SpaceDim);

  //set udotdelrho to zero.  the fortran is additive
  a_uDotDelRho.setVal(0.);
  BaseFab<Real>& regUDotDelRho = a_uDotDelRho.getSingleValuedFAB();
  //compute udotdelu everywhere as regular.
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceRho = a_faceRho[faceDir];
      const EBFaceFAB& faceVel = a_faceVel[faceDir];

      const BaseFab<Real>& regFaceRho = faceRho.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceVel.getSingleValuedFAB();

      //this does udotdelvel += 0.5*(uhigh+ulow)*(velhigh-vellow)/dx (non-conservative)
      //this does udotdelrho +=   (velhigh*rhohigh-vellow*rholow)/dx (    conservative)
      //where high == ivec + half*e^facedir
      FORT_ADVECTIVEF(CHF_FRA(regUDotDelRho),
                      CHF_CONST_FRA(regFaceRho),
                      CHF_CONST_FRA1(regFaceVel, 0),
                      CHF_CONST_INT(faceDir),
                      CHF_CONST_INT(ncomp),
                      CHF_CONST_REAL(m_dx[faceDir]),
                      CHF_BOX(a_box),
                      CHF_INT(s_doingVel));
    }

  //update the irregular vofs
  for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      //the set can be bigger than the box for performance reasons.
      if (a_box.contains(vof.gridIndex()))
        {
          // if (s_verbose && vof.gridIndex()[0]==12 && vof.gridIndex()[1]==8){pout() << "nonconservative advection " << vof << endl;}
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              //udelrho was set in regular update.  we reset it
              // to zero and recalc.
              Real uDelRhoPt = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, idir, Side::Lo);
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, idir, Side::Hi);
                  Real rhoLo = 0.0;
                  Real velLo = 0.0;
                  if (facesLo.size() > 0)
                    {
                      for (int iface = 0; iface < facesLo.size(); iface++)
                        {
                          rhoLo += a_faceRho[idir](facesLo[iface], ivar);
                          velLo += a_faceVel[idir](facesLo[iface], 0);
                        }
                      velLo /=  facesLo.size();
                      rhoLo /=  facesLo.size();
                      ch_flops()+=facesLo.size()*2+2;
                    }
                  else
                    {
                      const BaseIVFAB<Real>& coveredRhoLo = *a_coveredRhoLo[idir];
                      const BaseIVFAB<Real>& coveredVelLo = *a_coveredVelLo[idir];
                      rhoLo = coveredRhoLo(vof, ivar);
                      velLo = coveredVelLo(vof, 0);
                    }
                  Real rhoHi = 0.0;
                  Real velHi = 0.0;
                  if (facesHi.size() > 0)
                    {
                      for (int iface = 0; iface < facesHi.size(); iface++)
                        {
                          rhoHi += a_faceRho[idir](facesHi[iface], ivar);
                          velHi += a_faceVel[idir](facesHi[iface], 0);
                        }
                      velHi /=  facesHi.size();
                      rhoHi /=  facesHi.size();
                      ch_flops()+=facesHi.size()*2+2;
                    }
                  else
                    {
                      const BaseIVFAB<Real>& coveredRhoHi = *a_coveredRhoHi[idir];
                      const BaseIVFAB<Real>& coveredVelHi = *a_coveredVelHi[idir];
                      rhoHi = coveredRhoHi(vof, ivar);
                      velHi = coveredVelHi(vof, 0);
                    }

                  Real velAve = 0.5*(velHi + velLo);
                  Real rhoDiff = rhoHi - rhoLo;
                  uDelRhoPt += velAve*rhoDiff/m_dx[idir];
                } //end loop over directions
              a_uDotDelRho(vof, ivar) = uDelRhoPt;
              // if (s_verbose && vof.gridIndex()[0]==12 && vof.gridIndex()[1]==8){pout() << "   uDotDelRho " << uDelRhoPt << endl;}
            }//end loop over variables
            ch_flops()+=ncomp*SpaceDim*6;
        } //if vof is in this box
    }//end loop over vofs
}
/****************/
void
EBAdvectPatchIntegrator::
advectiveDerivative(EBCellFAB&                       a_uDotDelRho,
                    const EBFluxFAB&                 a_faceRho,
                    const EBFluxFAB&                 a_faceVel,
                    const Box&                       a_box)
{
  CH_TIME("EBPatchAdvect::advectiveDerivative");
  int ncomp = a_faceRho.nComp();


  //set udotdelrho to zero.  the fortran is additive
  a_uDotDelRho.setVal(0.);
  BaseFab<Real>& regUDotDelRho = a_uDotDelRho.getSingleValuedFAB();
  //compute udotdelu everywhere as regular.
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceRho = a_faceRho[faceDir];
      const EBFaceFAB& faceVel = a_faceVel[faceDir];

      const BaseFab<Real>& regFaceRho = faceRho.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceVel.getSingleValuedFAB();

      //this does udotdelvel += 0.5*(uhigh+ulow)*(velhigh-vellow)/dx (non-conservative)
      //this does udotdelrho +=   (velhigh*rhohigh-vellow*rholow)/dx (    conservative)
      //where high == ivec + half*e^facedir
      FORT_ADVECTIVEF(CHF_FRA(regUDotDelRho),
                      CHF_CONST_FRA(regFaceRho),
                      CHF_CONST_FRA1(regFaceVel, 0),
                      CHF_CONST_INT(faceDir),
                      CHF_CONST_INT(ncomp),
                      CHF_CONST_REAL(m_dx[faceDir]),
                      CHF_BOX(a_box),
                      CHF_INT(s_doingVel));
    }

}

#include "NamespaceFooter.H"
