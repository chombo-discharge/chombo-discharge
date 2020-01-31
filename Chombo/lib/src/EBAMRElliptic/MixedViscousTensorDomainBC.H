#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MIXEDVISCOUSTENSORDOMAINBC_H_
#define _MIXEDVISCOUSTENSORDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCFuncEval.H"
#include "ViscousBaseDomainBC.H"
#include "NamespaceHeader.H"

///
/**
 */
class MixedViscousTensorDomainBC: public ViscousBaseDomainBC
{
public:
  ///
  /**
   */
  MixedViscousTensorDomainBC()
  {
    ///set default start length to huge so it defaults to Neumann
    m_startDirichlet = 1.23456789e10;
  }

  ///
  /**
   */
  virtual ~MixedViscousTensorDomainBC()
  {;}


  ///
  /**
   */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);
  
  
  void setStartDirichlet(const Real& a_startDirichlet)
  {
    m_startDirichlet = a_startDirichlet;
  }

  virtual void 
  fillVelGhost(FArrayBox&     a_phi,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous)
  {

    Box grownBox = a_valid;
    grownBox.grow(1);
    for (int idir=0; idir<CH_SPACEDIM; ++idir)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            Box choppedBox = grownBox;
            choppedBox.grow(idir,-1);
            Box toRegion = adjCellBox(choppedBox, idir, sit(), 1);

            for (BoxIterator bit(toRegion); bit.ok(); ++bit)
              {
                const IntVect& iv = bit();
                //fake vof just to get the location
                VolIndex vof(iv, 0);
                RealVect loc = EBArith::getVoFLocation(vof, a_dx, RealVect::Zero);
                int isign = sign(sit());
                IntVect ivneigh = iv - isign*BASISV(idir);
                RealVect value = bcvaluefunc(loc, idir, sit());
                if(a_homogeneous) value = RealVect::Zero;
                for(int comp = 0; comp < SpaceDim; comp++)
                  {
                    if(loc[0] < m_startDirichlet)
                      {
                        //neum
                        a_phi(iv, comp) = a_phi(ivneigh, comp)  + a_dx*value[comp];
                      }
                    else
                      {
                        //diri
                        a_phi(iv, comp) = -a_phi(ivneigh, comp)  + 2.*value[comp];
                      }
                  }
              }
          } 
      }//end loop over directions
  }

private:

  Real m_startDirichlet;

};

#include "NamespaceFooter.H"
#endif
