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
#include "MOLPhysics.H"
#include "LoHiCenter.H"
#include "GodunovUtilitiesF_F.H"

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

/// Copy our values into one of us
void MOLPhysics::copyTo(MOLPhysics* a_mol) const
{
  CH_assert(a_mol != NULL);
  CH_assert(m_isBCSet);

  a_mol->setPhysIBC(m_bc);
  // added by petermc, 8 Jul 2009
  if (fourthOrderArtificialViscosityIsDefined())
      a_mol->setFourthOrderArtificialViscosityParameter(
          getFourthOrderArtificialViscosityParameter());
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
  CH_TIME("MOLPhysics::artVisc");
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

void MOLPhysics::soundSpeed(FArrayBox& a_speed,
                            const FArrayBox& a_U,
                            const Box&       a_box)
{
  MayDay::Abort("MOLPhysics::soundSpeed function not implemented");
}


Interval MOLPhysics::vectorFluxInterval() const
{
  // Default is the empty interval.
  return Interval();
}

#include "NamespaceFooter.H"
