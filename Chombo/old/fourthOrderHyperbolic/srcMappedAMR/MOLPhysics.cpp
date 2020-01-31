#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Box.H"
#include "CHArray.H"
#include "MOLPhysics.H"
#include "LoHiCenter.H"
// #include "GodunovPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "PhysMappedIBC.H"
#include "MOLPhysicsMappedArtViscF_F.H"

#include "NamespaceHeader.H"

MOLPhysics::MOLPhysics()
{
  m_isDefined = false;
  m_isBCSet = false;
  m_bc = NULL;
  m_useFourthOrderArtificialViscosity = false;
}

PhysIBC* MOLPhysics::getPhysIBC() const
{
  CH_assert(m_isBCSet);
  return m_bc;
}

void MOLPhysics::setPhysIBC(PhysIBC* a_bc)
{
  // Delete old boundary condition object - if any
  if (m_bc != NULL)
  {
    delete m_bc;
  }

  // Store new boundary condition object
  m_bc = a_bc->new_physIBC();

  // just in case we're re-defining the BC
  if (m_isDefined)
  {
    m_bc->define(m_domain, m_dx);
  }

  m_isBCSet = true;
}

MOLPhysics::~MOLPhysics()
{
  if (m_bc != NULL)
  {
    delete m_bc;
  }
}

void MOLPhysics::define(const ProblemDomain& a_domain,
                        const Real&          a_dx)
{
  m_domain    = a_domain;
  m_dx        = a_dx;
  m_isDefined = true;

  if (m_bc != NULL)
  {
    m_bc->define(m_domain, m_dx);
  }

  m_util.define(m_domain, m_dx);
}

void MOLPhysics::setCurrentBox(const Box& a_currentBox)
{
  // Do nothing but assert the object has been defined
  CH_assert(isDefined());
}
void MOLPhysics::setCurrentCoordSys(const NewCoordSys* a_coordSys)
{
  // Do nothing but assert the object has been defined
  CH_assert(isDefined());
}
void MOLPhysics::setCurrentTime(const Real& a_currentTime)
{
  // Do nothing but assert the object has been defined
  CH_assert(isDefined());
}
void MOLPhysics::setFourthOrderArtificialViscosityParameter(const Real& a_M0sq)
{
  m_M0sq = a_M0sq;
  m_useFourthOrderArtificialViscosity = true;
}
bool MOLPhysics::fourthOrderArtificialViscosityIsDefined() const
{
  return m_useFourthOrderArtificialViscosity;
}
Real MOLPhysics::getFourthOrderArtificialViscosityParameter() const
{
  return m_M0sq;
}
int MOLPhysics::densityIndex()
{
  return -1;
}

void MOLPhysics::getFlux(FArrayBox&       a_flux,
                         const FArrayBox& a_WHalf,
                         const int&       a_dir,
                         const Box&       a_box)
{
  MayDay::Error("MOLPhysics::getFlux:  Default implementation called - this should never happen");
}

// Set face-averaged primitive state on boundary faces
void MOLPhysics::primBC(FArrayBox&             a_WGdnv,
                        const FArrayBox&       a_WLeft,
                        const FArrayBox&       a_WRight,
                        const FArrayBox&       a_W,
                        const FArrayBox *const a_unitNormalBasisPtr,
                        const Real&            a_time,
                        const int&             a_dir)
{
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
  static_cast<PhysMappedIBC*>(m_bc)->primBC(
    a_WGdnv, shiftWLeft , a_W, a_unitNormalBasisPtr, velocityInterval(), a_dir,
    Side::Hi, a_time);
  static_cast<PhysMappedIBC*>(m_bc)->primBC(
    a_WGdnv, shiftWRight, a_W, a_unitNormalBasisPtr, velocityInterval(), a_dir,
    Side::Lo, a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

bool MOLPhysics::isDefined() const
{
  return m_isDefined;
}

void MOLPhysics::artVisc(FArrayBox&       a_F,
                         const FArrayBox& a_U,
                         const Real&      a_artificialViscosity,
                         const Real&      a_currentTime,
                         const int&       a_dir,
                         const Box&       a_box)
{
  CH_assert(a_U.box().contains(a_box)); // a_box:  valid cells

  // Set faceBox0 to all a_dir-faces of a_box.
  Box faceBox0 = a_box;
  faceBox0.surroundingNodes(a_dir);

  CH_assert(a_F.box().contains(faceBox0));

  // In 4th-order case, bx1 will be shrunken by 1 in each transverse direction
  // after we find 3-point minimum in that direction.
  Box bx1 = grow(a_box, 1);

  // Need primitive variables on bx1inDomain.
  Box bx1inDomain(bx1);
  bx1inDomain &= m_domain;

  // Get the primitive variables from the conserved variables (as needed).
  int numPrim = numPrimitives();
  FArrayBox W(bx1inDomain, numPrim);
  /*
    Find the primitive variables
    W on bx1inDomain = grow(a_box, 1) & m_domain
    using a_U on bx1inDomain = grow(a_box, 1) & m_domain.
   */
  consToPrim(W, a_U, bx1inDomain);

  /*
    Compute the divergence of the velocity
    divu on faceBox0 = all a_dir-faces of a_box
    using W on bx1inDomain = grow(a_box, 1) & m_domain.
  */
  FArrayBox divu(faceBox0, 1);
  Interval velInt = velocityInterval();
  m_util.divVel(divu, W, velInt, a_dir, faceBox0);
  // If using fourth-order artificial viscosity, apply the nonlinear operator to
  // divu.
  if (m_useFourthOrderArtificialViscosity)
    {
      // m_util.divVelHO(divu, W, a_dir, faceBox0, this);
      CH_assert(fourthOrderArtificialViscosityIsDefined());

      // Compute cell-centered (bulk modulus)/(density).
      int bulkIndex = bulkModulusIndex();
      int densIndex = densityIndex();
      Real M0sq = getFourthOrderArtificialViscosityParameter();
      FArrayBox csq1(bx1inDomain, 1);
      FArrayBox csq2(bx1inDomain, 1);
      /*
        Set csq1 = W[bulkIndex] / W[densIndex]
        on all cells of bx1inDomain = grow(a_box, 1) & m_domain.
       */
      csq1.setVal(0.);
      csq2.setVal(0.);
      csq1.copy(W, bx1inDomain, bulkIndex, bx1inDomain, 0, 1);
      csq1.divide(W, bx1inDomain, bx1inDomain, densIndex, 0, 1);
      Box hiBox,loBox,centerBox,entireBox;
      int hasLo,hasHi;
      FArrayBox* csqin_ptr = &csq1;
      FArrayBox* csqout_ptr = &csq2;

      // Compute min of csq in transverse direction(s).
      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          if (dir != a_dir)
            {
              FArrayBox& csqin = *csqin_ptr;
              FArrayBox& csqout = *csqout_ptr;
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         bx1, m_domain,dir);

              FORT_MIN3PTSF(CHF_FRA1(csqout,0),
                            CHF_CONST_FRA1(csqin,0),
                            CHF_CONST_INT(dir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));
              bx1.grow(dir, -1);
              FArrayBox* csqtmp = csqin_ptr;
              csqin_ptr = csqout_ptr;
              csqout_ptr = csqtmp;
            }
        }
      // bx1dir = valid cells + 1 ghost layer in direction a_dir
      Box bx1dir(a_box);
      bx1dir.grow(a_dir, 1);
      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     bx1dir, m_domain, a_dir);

      FArrayBox& csq = *csqin_ptr;
      /*
        Modify divu on all a_dir-faces of a_box
        using csq on all cells of bx1inDomain = grow(a_box, 1) & m_domain.
       */
      FORT_HIGHORDERDIVCO(CHF_FRA1(divu,0),
                          CHF_CONST_FRA1(csq,0),
                          CHF_CONST_INT(a_dir),
                          CHF_CONST_REAL(M0sq),
                          CHF_BOX(loBox),
                          CHF_CONST_INT(hasLo),
                          CHF_BOX(hiBox),
                          CHF_CONST_INT(hasHi),
                          CHF_BOX(centerBox));
    }

  // Change fluxes due to artificial viscosity on
  // faceInteriorBox = all a_dir-faces of a_box that are not on boundaries.
  /*
    Modify a_F
    on faceInteriorBox = all a_dir-faces of a_box that are not on boundaries
    using a_U on all cells of grow(a_box, BASISV(a_dir)) & m_domain
    and divu on faceInteriorBox.
   */
  Box faceInteriorBox(a_box);
  faceInteriorBox.grow(a_dir, 1);
  faceInteriorBox &= m_domain;
  faceInteriorBox.grow(a_dir, -1);
  faceInteriorBox.surroundingNodes(a_dir);
  m_util.artificialViscosity(a_F, a_U,
                             divu, a_artificialViscosity, a_dir,
                             faceInteriorBox);

  // Change fluxes due to artificial viscosity on the boundary faces
  m_bc->artViscBC(a_F, a_U, divu, a_dir, a_currentTime);
}

void MOLPhysics::mappedArtVisc(FluxBox&         a_NtF,
                               const FArrayBox& a_U,
                               const FluxBox&   a_N,
                               const FArrayBox& a_J,
                               const FluxBox&   a_unitNormalBasis,
                               const Real&      a_alpha,
                               const Real&      a_currentTime,
                               const Box&       a_box)
{
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

      static_cast<PhysMappedIBC*>(m_bc)->artViscBC(
        a_NtF[dir],
        Nctg,
        a_U,
        a_unitNormalBasis[dir],
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

void MOLPhysics::soundSpeed(FArrayBox& a_speed,
                            const FArrayBox& a_U,
                            const Box&       a_box)
{
  MayDay::Abort("MOLPhysics::soundSpeed function not implemented");
}

#include "NamespaceFooter.H"
