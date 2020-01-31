#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>

#include "EBInflowOutflowIBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "ParmParse.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "EBInflowOutflowF_F.H"
#include "EBPatchPolytropicF_F.H"
/****************************/
/****************************/
EBInflowOutflowIBC::
EBInflowOutflowIBC(const Real&     a_gamma,
                   const Real&     a_cv,
                   const Real&     a_tempinflow,
                   const Real&     a_tempinterior,
                   const Real&     a_presinflow,
                   const Real&     a_presinterior,
                   const Real&     a_machinflow,
                   const int&      a_flowdir)        
  :EBPhysIBC()
{
  m_gamma         = a_gamma        ;
  m_cv            = a_cv           ;
  m_tempinflow    = a_tempinflow   ;
  m_tempinterior  = a_tempinterior ;
  m_presinflow    = a_presinflow   ;
  m_presinterior  = a_presinterior ;
  m_machinflow    = a_machinflow   ;
  m_flowdir       = a_flowdir      ;

  //sets fortran commons
  FORT_SETINFLOWOUTFLOW(CHF_CONST_REAL(a_cv),
                        CHF_CONST_REAL(a_gamma),
                        CHF_CONST_REAL(a_tempinflow),
                        CHF_CONST_REAL(a_tempinterior),
                        CHF_CONST_REAL(a_machinflow),
                        CHF_CONST_REAL(a_presinflow),
                        CHF_CONST_REAL(a_presinterior));

  m_isFortranCommonSet = true;
  m_isDefined = false;
}
/****************************/
/****************************/
void
EBInflowOutflowIBC::
fluxBC(EBFluxFAB&            a_flux,
       const EBCellFAB&      a_primCenter,
       const EBCellFAB&      a_primExtrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_dir)
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);
 CH_assert(!m_domain.isPeriodic(a_dir));

  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;


  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      bool atInflow  = ((a_dir == m_flowdir) && (a_side == Side::Lo));
      bool atOutflow = ((a_dir == m_flowdir) && (a_side == Side::Hi));

      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box boundaryBox = bdryBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-isign);
      BaseFab<Real>& regFlux = a_flux[a_dir].getSingleValuedFAB();
      const BaseFab<Real>& regPrimExtrap = a_primExtrap.getSingleValuedFAB(); // 

      Box region = a_primExtrap.box();
      region &= a_ebisBox.getDomain();

      regFlux.shiftHalf(a_dir,-isign);

      // Set the boundary fluxes
      /**/
      if(atInflow)
        {
          FORT_EBIOINFLOWBC( CHF_FRA(regFlux),
                             CHF_CONST_FRA(regPrimExtrap),
                             CHF_CONST_INT(isign),
                             CHF_CONST_INT(a_dir),
                             CHF_BOX(boundaryBox));
        } 
      else if(atOutflow)
        {
          FORT_EBIOOUTFLOWBC(CHF_FRA(regFlux),
                             CHF_CONST_FRA(regPrimExtrap),
                             CHF_CONST_INT(isign),
                             CHF_CONST_INT(a_dir),
                             CHF_BOX(boundaryBox));
        }
      else
        {
          FORT_EBIOSLIPWALLBC(CHF_FRA(regFlux),
                              CHF_CONST_FRA(regPrimExtrap),
                              CHF_CONST_INT(isign),
                              CHF_CONST_INT(a_dir),
                              CHF_BOX(boundaryBox));

        }

      // Shift returned fluxes to be face centered
      regFlux.shiftHalf(a_dir,isign);
      //now for the multivalued cells.  Since it is pointwise,
      //the regular calc is correct for all single-valued cells.
      IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
      for(VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Vector<Real> qgdnv(QNUM);
          Vector<Real> fluxv(FNUM);
          for(int ivar = 0; ivar < QNUM; ivar++)
            {
              qgdnv[ivar] = a_primExtrap(vof, ivar);
            }
          Vector<FaceIndex> bndryFaces =
            a_ebisBox.getFaces(vof, a_dir, a_side);
          for(int iface= 0; iface < bndryFaces.size(); iface++)
            {
              /**/
              if(atInflow)
                {
                  FORT_POINTEBIOINFLOWBC(CHF_VR(fluxv),
                                          CHF_VR(qgdnv),
                                          CHF_CONST_INT(a_dir),
                                          CHF_CONST_INT(isign)
                                          );
                } 
              else if(atOutflow)
                {
                  FORT_POINTEBIOOUTFLOWBC(CHF_VR(fluxv),
                                          CHF_VR(qgdnv),
                                          CHF_CONST_INT(a_dir),
                                          CHF_CONST_INT(isign)
                                          );
                }
              else
                {
                  FORT_POINTEBIOSLIPWALLBC(CHF_VR(fluxv),
                                           CHF_VR(qgdnv),
                                           CHF_CONST_INT(a_dir),
                                           CHF_CONST_INT(isign)
                                           );
                }

              const FaceIndex& face = bndryFaces[iface];
              for(int ivar = 0; ivar < FNUM; ivar++)
                {
                  a_flux[a_dir](face, ivar) = fluxv[ivar];
                }
            }
        }
    }
}
/****************************/
/****************************/
void
EBInflowOutflowIBC::define(const ProblemDomain&  a_domain,
                         const RealVect& a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx[0];
 CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);
  m_isDefined = true;
}
/****************************/
/****************************/
void
EBInflowOutflowIBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

  // Iterator of all grids in this level
  for(DataIterator dit = a_conState.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      // Storage for current grid
      EBCellFAB& conFAB = a_conState[dit()];

      BaseFab<Real>& regConFAB = conFAB.getSingleValuedFAB();
      // Box of current grid
      Box uBox = regConFAB.box();
      uBox &= m_domain;
      // Set up initial condition in this grid
      /**/
      FORT_EBIOINFLOWOUTFLOWINIT(CHF_CONST_FRA(regConFAB),
                                 CHF_CONST_REAL(m_dx),
                                 CHF_BOX(uBox));
      /**/

      //now for the multivalued cells.  Since it is pointwise,
      //the regular calc is correct for all single-valued cells.
      IntVectSet ivs = ebisBox.getMultiCells(uBox);
      for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          RealVect momentum;
          Real energy, density;
          /**/
          FORT_POINTEBIOINFLOWOUTFLOWINIT(CHF_REAL(density),
                                          CHF_REALVECT(momentum),
                                          CHF_REAL(energy),
                                          CHF_CONST_INTVECT(iv),
                                          CHF_CONST_REAL(m_dx));
          /**/
          conFAB(vof, CRHO) = density;
          conFAB(vof, CENG) = energy;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              conFAB(vof, CMOMX+idir) = momentum[idir];
            }
        }//end loop over multivalued cells
    } //end loop over boxes
}
/****************************/
/****************************/
EBInflowOutflowIBC::
~EBInflowOutflowIBC()
{
}
/****************************/
/****************************/
