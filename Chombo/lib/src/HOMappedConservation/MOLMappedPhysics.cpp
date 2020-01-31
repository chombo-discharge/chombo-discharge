#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MOLMappedPhysics.H"
#include "Box.H"
#include "CHArray.H"
#include "LoHiCenter.H"
// #include "GodunovPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "PhysMappedIBC.H"
#include "MOLPhysicsMappedArtViscF_F.H"

#include "NamespaceHeader.H"

MOLMappedPhysics::MOLMappedPhysics() : 
  MOLPhysics(), 
  m_time(0.0), 
  m_coordSys(NULL),
  m_haveUnitNormalBasis(false)
{
}


MOLMappedPhysics::~MOLMappedPhysics()
{
}


// Given a MOLMappedPhysics pointer (thus could be a derived class), 
// copy over the member data only we know about.
void MOLMappedPhysics::copyTo(MOLMappedPhysics* a_newMOLMP) const
{
  CH_assert(a_newMOLMP != NULL);
  MOLPhysics::copyTo(a_newMOLMP);

  // Copy our derived class members over
  if (m_coordSys != NULL)
    a_newMOLMP->setCurrentCoordSys(m_coordSys);
  a_newMOLMP->setCurrentTime(m_time);
}


void MOLMappedPhysics::setCurrentTime(const Real& a_time)
{
  m_time = a_time;
}


void MOLMappedPhysics::setCurrentCoordSys(const NewCoordSys* a_coordSys)
{
  m_coordSys = const_cast<NewFourthOrderCoordSys*>(
      dynamic_cast<const NewFourthOrderCoordSys*>(a_coordSys));
  CH_assert(m_coordSys != NULL);
}


// Set face-averaged primitive state on boundary faces
void MOLMappedPhysics::primBC(FArrayBox&             a_WGdnv,
                              const FArrayBox&       a_WLeft, 
                              const FArrayBox&       a_WRight, 
                              const FArrayBox&       a_W, 
                              const Real&            a_time, 
                              const int&             a_dir)
{
  CH_TIME("MOLMappedPhysics::primBC");
  CH_assert(isDefined());

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  // Set the primitive state on boundary faces
  CH_assert(m_haveUnitNormalBasis);
  const FArrayBox* const unitNormalDirFabPtr = &((*m_unitNormalBasisPtr)[a_dir]);
  static_cast<PhysMappedIBC*>(m_bc)->primBC(
    a_WGdnv, shiftWLeft , a_W, unitNormalDirFabPtr, velocityInterval(), a_dir,
    Side::Hi, a_time);
  static_cast<PhysMappedIBC*>(m_bc)->primBC(
    a_WGdnv, shiftWRight, a_W, unitNormalDirFabPtr, velocityInterval(), a_dir,
    Side::Lo, a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}


void MOLMappedPhysics::setUnitNormalBasisPtr(const FluxBox* const a_unitNormalBasisPtr)
{
  m_unitNormalBasisPtr = (FluxBox*) a_unitNormalBasisPtr;
  m_haveUnitNormalBasis = true;
}

void MOLMappedPhysics::mappedArtVisc(FluxBox&         a_NtF,
                                     const FArrayBox& a_U,
                                     const FluxBox&   a_N, 
                                     const FArrayBox& a_J, 
                                     const Real&      a_alpha, 
                                     const Real&      a_currentTime, 
                                     const Box&       a_box)
{
  CH_TIME("MOLMappedPhysics::mappedArtVisc");
  CH_assert(fourthOrderArtificialViscosityIsDefined());
  CH_assert(a_NtF.box().contains(a_box));

  Box bx1 = grow(a_box, 1);
  Box bx1inDomain(bx1);
  bx1inDomain &= m_domain;
  CH_assert(a_U.box().contains(bx1inDomain));

  // Cell-centered boxes providing marking low and high side of the domain
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  // Get the primitive variables from the conserved variables (as needed).
  int numPrim = numPrimitives();
  FArrayBox W(bx1inDomain, numPrim);
  consToPrim(W, a_U, bx1inDomain);
  FArrayBox vel(velocityInterval(), W);  // Alias

//--Precompute gradients of vel & U in all cells.  Note that the dyads have
//--numComp*SpaceDim for all gradient directions.  gradVel is needed for
//--divergence of the velocity and gradU for the artificial viscosity.

  const int numVelComp = SpaceDim;
  // Components stored as (velocity direction, gradient direction) in Fortran
  // ordering
  const int velCompStride = 1;
  CHArray<Real, SpaceDim+1, ArRangeCol> gradVel(numVelComp*SpaceDim,
                                                bx1inDomain);
  gradVel = -1.;  // Not strictly required but avoids computation with
                  // uninitialized values

  const int numUComp = numConserved();
  // Components stored as (gradient direction, conserved variable) in Fortran
  // ordering
  const int UCompStride = SpaceDim;
  CHArray<Real, SpaceDim+1, ArRangeCol> gradU(SpaceDim*numConserved(),
                                              bx1inDomain);
  gradU = -1.;    // Not strictly required but avoids computation with
                  // uninitialized values

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // LHCbox needs to be grown by 1 in dir and not intersected with the
      // domain to satisfy the requirements of loHiCenter.
      // In the other directions, it needs to be grown by 1 and intersected with
      // the domain.  This is so we have a centered understanding of gradients
      // tangential to a face at the edge of a box (but not the edge of the
      // domain).  Because we only need tangential gradients outside a_box, the
      // velocity and U are still only required within 1 ghost cell of a_box.
      Box LHCbox(bx1);
      LHCbox.grow(dir, -1);
      LHCbox &= m_domain;
      LHCbox.grow(dir, 1);
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, LHCbox,
                 m_domain, dir);

      const int velCompBegin = dir*numVelComp;
      FORT_CELLGRADDIR(
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradVel),
        CHF_CONST_FRA(vel),
        CHF_BOX(loBox),
        CHF_CONST_INT(hasLo),
        CHF_BOX(hiBox),
        CHF_CONST_INT(hasHi),
        CHF_BOX(centerBox),
        CHF_CONST_INT(dir),
        CHF_CONST_INT(numVelComp),
        CHF_CONST_INT(velCompBegin),
        CHF_CONST_INT(velCompStride),
        CHF_CONST_REAL(m_dx));

      const int UCompBegin = dir;
      FORT_CELLGRADDIR(
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradU),
        CHF_CONST_FRA(a_U),
        CHF_BOX(loBox),
        CHF_CONST_INT(hasLo),
        CHF_BOX(hiBox),
        CHF_CONST_INT(hasHi),
        CHF_BOX(centerBox),
        CHF_CONST_INT(dir),
        CHF_CONST_INT(numUComp),
        CHF_CONST_INT(UCompBegin),
        CHF_CONST_INT(UCompStride),
        CHF_CONST_REAL(m_dx));
    }

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box faceBox0 = a_box;
      faceBox0.surroundingNodes(dir);

//--Store a contiguous N

      const int numNComp = SpaceDim*SpaceDim;
      const int zeroVal = 0;
      CHArray<Real, SpaceDim+1, ArRangeCol> Nctg(numNComp, faceBox0);
      FORT_REVERSEFABCOMPONENTSTRIDE(CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                                     CHF_CONST_FRA(a_N[dir]),
                                     CHF_BOX(faceBox0),
                                     CHF_CONST_INT(zeroVal),
                                     CHF_CONST_INT(zeroVal),
                                     CHF_CONST_INT(numNComp));

//--Compute cell-centered (bulk modulus)/(density), i.e., c^2/gamma

      FArrayBox csqA(bx1inDomain, 1);
      FArrayBox csqB(bx1inDomain, 1);

      csqA.setVal(0.);
      csqA.setVal(0.);
      csqA.copy(W, bx1inDomain, bulkModulusIndex(), bx1inDomain, 0, 1);
      csqA.divide(W, bx1inDomain, densityIndex(), 0, 1);

      FArrayBox* csqIptr = &csqA;
      FArrayBox* csqOptr = &csqB;

      // Compute min of csq in transverse direction(s).
      Box csqBox(bx1);
      for (int trDir = 0; trDir != SpaceDim; ++trDir)
        {
          if (trDir != dir)
            {
              FArrayBox& csqI = *csqIptr;
              FArrayBox& csqO = *csqOptr;
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         csqBox, m_domain, trDir);

              FORT_MIN3PTSF(CHF_FRA1(csqO, 0),
                            CHF_CONST_FRA1(csqI, 0),
                            CHF_CONST_INT(trDir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));
              csqBox.grow(trDir, -1);
              FArrayBox* csqTmp = csqIptr;
              csqIptr = csqOptr;
              csqOptr = csqTmp;
            }
        }
      FArrayBox& csq = *csqIptr;

//--Face boxes for this direction

      Box bx1dir(a_box);
      bx1dir.grow(dir, 1);
      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, bx1dir,
                     m_domain, dir);

//--Compute the divergence of velocity

      FArrayBox divVel(faceBox0, 1);
      FORT_MAPPEDDIVVEL(CHF_FRA1(divVel, 0),
                        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradVel),
                        CHF_CONST_FRA(vel),
                        CHF_CONST_FRA1(a_J, 0),
                        CHF_BOX(loBox),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(hiBox),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(centerBox),
                        CHF_CONST_INT(dir),
                        CHF_CONST_REAL(m_dx));

//--Compute the physical cell spacing across the faces

      FArrayBox dxFace(faceBox0, 1);
      FORT_PHYSICALCELLSPACINGONFACE(CHF_FRA1(dxFace, 0),
                                     CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                                     CHF_CONST_FRA1(a_J, 0),
                                     CHF_BOX(loBox),
                                     CHF_CONST_INT(hasLo),
                                     CHF_BOX(hiBox),
                                     CHF_CONST_INT(hasHi),
                                     CHF_BOX(centerBox),
                                     CHF_CONST_INT(dir),
                                     CHF_CONST_REAL(m_dx));

//--Compute the flux due to artificial viscosity on the faces of the cells in
//--*computational* space (i.e., they have already been multiplied by a row of
//--N^T).  Only interior faces are affected.

      // Need to set boundary values to zero in case they are not modified
      if (hasLo)
        {
          a_NtF[dir].setVal(0., loBox, 0, a_NtF.nComp());
        }
      if (hasHi)
        {
          a_NtF[dir].setVal(0., hiBox, 0, a_NtF.nComp());
        }

      const int hasLoHiFalse = 0;
      const Real beta = getFourthOrderArtificialViscosityParameter();
      FORT_MAPPEDARTVISC(CHF_FRA(a_NtF[dir]),
                         CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                         CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradU),
                         CHF_CONST_FRA(a_U),
                         CHF_CONST_FRA1(divVel, 0),
                         CHF_CONST_FRA1(csq, 0),
                         CHF_CONST_FRA1(a_J, 0),
                         CHF_CONST_FRA1(dxFace, 0),
                         CHF_CONST_REAL(a_alpha),
                         CHF_CONST_REAL(beta),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLoHiFalse),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasLoHiFalse),
                         CHF_BOX(centerBox),
                         CHF_CONST_INT(dir),
                         CHF_CONST_REAL(m_dx));

//--Change fluxes due to artificial viscosity on the boundary faces

      CH_assert(m_haveUnitNormalBasis);
      static_cast<PhysMappedIBC*>(m_bc)->artViscBC(
        a_NtF[dir],
        Nctg,
        a_U,
        (*m_unitNormalBasisPtr)[dir],
        divVel,
        csq,
        dxFace,
        velocityInterval(),  // Expected to be the same as momentum
        a_alpha,
        beta,
        loBox,
        hasLo,
        hiBox,
        hasHi,
        dir);
    }
}

void MOLMappedPhysics::getAux(FArrayBox& a_auxFAb,
                              const FArrayBox& a_UFab)
{
  // default: do nothing
}

Vector<string> MOLMappedPhysics::auxNames()
{
  // default: no aux data, so return empty Vector
  return Vector<string>();
}

#include "NamespaceFooter.H"
