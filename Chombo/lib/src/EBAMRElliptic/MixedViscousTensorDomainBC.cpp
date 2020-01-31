#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"
#include "EBArith.H"

#include "MixedViscousTensorDomainBC.H"
#include "TensorCFInterp.H"
#include "ViscousTensorOpF_F.H"
#include "EBViscousTensorOp.H"
#include "EBViscousTensorOpF_F.H"
#include "DirichletViscousTensorDomainBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "NamespaceHeader.H"



void
MixedViscousTensorDomainBC::
getFaceFlux(BaseFab<Real>&        a_flux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  Box box = a_flux.box();
  BaseFab<Real> diriflux(box, a_flux.nComp());
  BaseFab<Real> neumflux(box, a_flux.nComp());
  DirichletViscousTensorDomainBC diriBC;
  NeumannViscousTensorDomainBC   neumBC;
  diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
  neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda); 
  diriBC.setValue(0.);
  neumBC.setValue(0.);
  diriBC.getFaceFlux(diriflux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
  neumBC.getFaceFlux(neumflux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);

  for (int comp=0; comp<a_phi.nComp(); comp++)
    {

      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          for(int comp = 0; comp < a_flux.nComp(); comp++)
            {
              RealVect point = EBArith::getIVLocation(iv,a_dx,a_probLo);
              Real fluxval;
              if(point[0] > m_startDirichlet)
                {
                  fluxval = diriflux(iv, comp);
                }
              else
                {
                  fluxval = neumflux(iv, comp);
                }

              a_flux(iv,comp) = fluxval;
            }
        }
    }
}
/****/
void
MixedViscousTensorDomainBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  RealVect  loc = EBArith::getVofLocation(a_vof, a_dx, RealVect::Zero);
  if(loc[0] > m_startDirichlet)
    {
      DirichletViscousTensorDomainBC diriBC;
      diriBC.setValue(0.);
      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      NeumannViscousTensorDomainBC neumBC;
      neumBC.setValue(0.);
      neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}
#include "NamespaceFooter.H"
