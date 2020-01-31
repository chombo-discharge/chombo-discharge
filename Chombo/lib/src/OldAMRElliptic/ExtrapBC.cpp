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
#include "FArrayBox.H"
#include "REAL.H"
#include "SPACE.H"
#include "Tuple.H"
#include "Vector.H"
#include "LoHiSide.H"

#include "ExtrapBC.H"
#include "MayDay.H"
#include "ExtrapBCF_F.H"
#include "NamespaceHeader.H"

ExtrapBC::ExtrapBC() : BoxGhostBC()
{
}

ExtrapBC::~ExtrapBC()
{
}

void
ExtrapBC::fillBCValues(FArrayBox& a_neumfac,
                       FArrayBox& a_dircfac,
                       FArrayBox& a_inhmval,
                       Real a_dx,
                       const ProblemDomain& a_domain) const
{
  // this doesn't do anything because we don't need these arrays...
  // this function is just needed to complete the instantiation.
}


void
ExtrapBC::fillBCValues(FArrayBox& a_neumfac,
                       FArrayBox& a_dircfac,
                       FArrayBox& a_inhmval,
                       Real a_dx,
                       const Box& a_domain) const
{
  // this doesn't do anything because we don't need these arrays...
  // this function is just needed to complete the instantiation.
}


BoxGhostBC*
ExtrapBC::new_boxghostbc() const
{
  ExtrapBC* newop = new ExtrapBC();
  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in ExtrapBC::new_boxghostbc");
    }
  return static_cast<BoxGhostBC*>(newop);
}

// -------------------------------------------------------------
void
ExtrapBC::applyHomogeneousBCs(FArrayBox& a_state,
                              const ProblemDomain& a_domain,
                              Real a_dx) const
{
  applyExtrapBCs(a_state, a_domain, a_dx);
}

// -------------------------------------------------------------
void
ExtrapBC::applyHomogeneousBCs(FArrayBox& a_state,
                              const Box& a_domain,
                              Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyExtrapBCs(a_state, physdomain, a_dx);
}

// -------------------------------------------------------------
void
ExtrapBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                const ProblemDomain& a_domain,
                                Real a_dx) const
{
  applyExtrapBCs(a_state, a_domain, a_dx);
}

// -------------------------------------------------------------
void
ExtrapBC::applyInhomogeneousBCs(FArrayBox& a_state,
                                const Box& a_domain,
                                Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyExtrapBCs(a_state, physdomain, a_dx);
}



// -------------------------------------------------------------
void
ExtrapBC::applyExtrapBCs(FArrayBox& a_state, const Box& a_domain,
                         Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyExtrapBCs(a_state, a_domain, a_dx);
}


// -------------------------------------------------------------
void
ExtrapBC::applyExtrapBCs(FArrayBox& a_state,
                         const ProblemDomain& a_domain,
                         Real a_dx) const
{

  // don't do anything in periodic case
  if (!a_domain.isPeriodic(m_direction))
  {
    // first construct BC box
    Box bx= a_state.box();
    bx.grow(-1);
    Box bc_box;
    bool isbc;
    int idir = m_direction;
    if (m_side == Side::Lo)
      {
        isbc = (bx.smallEnd(idir) <= a_domain.domainBox().smallEnd(idir));
        if (isbc)
          {
            int ichop = a_domain.domainBox().smallEnd(m_direction);
            bc_box = a_state.box();
            bc_box.chop(m_direction, ichop);
          }
      }
    else if (m_side == Side::Hi)
      {
        isbc = (bx.bigEnd(idir) >= a_domain.domainBox().bigEnd(idir));
        if (isbc)
          {
            int ichop = a_domain.domainBox().bigEnd(m_direction)+1;
            Box chop_box = a_state.box();
            bc_box = chop_box.chop(m_direction,ichop);
          }
      }
    else
      {
        cerr << "DomainGhostBC::applyghostbc: bogus side" << endl;
        abort();
      }

    if (isbc)
      {
        int iSide = sign(m_side);
        int iDir = m_direction;

#ifndef NDEBUG
        CH_assert((m_side == Side::Hi)||(m_side == Side::Lo));
        Box biggerbox = bc_box;
        //state has to contain one cell on side of ghost cell
        if (m_side == Side::Lo)
          biggerbox.growHi(iDir, 1);
        else
          biggerbox.growLo(iDir, 1);


        CH_assert(a_state.box().contains(biggerbox));

        CH_assert(m_components.begin() >= 0);
        CH_assert(m_components.end() < a_state.nComp());
#endif
        int startcomp = m_components.begin();
        int endcomp = m_components.end();
        FORT_EXTRAPGHOSTBC(CHF_FRA(a_state),
                           CHF_BOX(bc_box),
                           CHF_CONST_INT(iDir),
                           CHF_CONST_INT(iSide),
                           CHF_CONST_REAL(a_dx),
                           CHF_CONST_INT(startcomp),
                           CHF_CONST_INT(endcomp));
      } // end if this is a bc

  }
}


ExtrapBC::ExtrapBC(int dir, Side::LoHiSide side, const Interval& a_comps)
    : BoxGhostBC(dir,side,a_comps)
{
}


ExtrapBC::ExtrapBC(int dir, Side::LoHiSide side)
  : BoxGhostBC(dir,side)
{
}



#include "NamespaceFooter.H"
