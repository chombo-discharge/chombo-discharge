/*!
  @file   jump_bcI.H
  @brief  Implementation of jump_bc.cpp
  @author Robert Marskar
  @date   Dec. 2017
*/

#ifndef _JUMP_BCI_
#define _JUMP_BCI_

#include "jump_bc.H"

#define DEBUG 0

inline
void jump_bc::compute_avg_jump(const BaseIVFAB<Real>& a_jump, const MFCellFAB& a_phi, const DataIndex& a_dit){

  //  const IntVectSet& ivs = m_avgJump[a_dit].get_ivs();

  const int comp   = 0;
  const int phase1 = 0;
  const int phase2 = 1;
  
  const BaseIVFAB<Real>& bco1 = m_bco[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& bco2 = m_bco[a_dit].get_ivfab(phase2);

  BaseIVFAB<Real>& avg1 = m_avgJump[a_dit].get_ivfab(phase1);
  BaseIVFAB<Real>& avg2 = m_avgJump[a_dit].get_ivfab(phase2);

  avg1.setVal(0.0);
  avg2.setVal(0.0);
  
  const EBISBox& ebisbox1 = a_phi.getPhase(phase1).getEBISBox();
  const EBISBox& ebisbox2 = a_phi.getPhase(phase2).getEBISBox();

  for (IVSIterator ivsit(m_ivs[a_dit]); ivsit.ok(); ++ivsit){
    const IntVect iv = ivsit();
    
    const Vector<VolIndex> vofs1 = ebisbox1.getVoFs(iv);
    const Vector<VolIndex> vofs2 = ebisbox2.getVoFs(iv);

    // Compute average jump condition

    Real totalJump = 0.0;

    if(ebisbox1.isMultiValued(iv)){
      Real totalArea = jump_bc::SAFETY;
      
      for (int v = 0; v < vofs1.size(); v++){
	const VolIndex vof = vofs1[v];

	const Real area = ebisbox1.bndryArea(vof);

	totalArea += area;
	totalJump += area*a_jump(vof,0);
      }
      totalJump = totalJump/totalArea;
    }
    else{
      totalJump = a_jump(vofs1[0], comp);
    }


    // Set the jump
    const VolIndex vof0(iv, comp);
    avg1(vof0, comp) = totalJump;
    avg2(vof0, comp) = totalJump;
  }
}

inline
void jump_bc::match_bc(BaseIVFAB<Real>&                  a_phibc,
		       const BaseIVFAB<Real>&            a_jump,
		       const MFCellFAB&                  a_phi,
		       const DataIndex&                  a_dit){
  const int comp   = 0;
  const int phase1 = 0;
  const int phase2 = 1;

  // First, compute the area-weighted jump coefficient
  jump_bc::compute_avg_jump(a_jump, a_phi, a_dit);
  
  const BaseIVFAB<Real>& avg_bco1    = m_avgBco[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& avg_bco2    = m_avgBco[a_dit].get_ivfab(phase2);
  
  const BaseIVFAB<Real>& avg_w1      = m_avgWeights[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& avg_w2      = m_avgWeights[a_dit].get_ivfab(phase2);
  
  const BaseIVFAB<VoFStencil>& avg_sten1 = m_avgStencils[a_dit].get_ivfab(phase1);
  const BaseIVFAB<VoFStencil>& avg_sten2 = m_avgStencils[a_dit].get_ivfab(phase2);

  const BaseIVFAB<Real>& avg_jump1   = m_avgJump[a_dit].get_ivfab(phase1);
  const BaseIVFAB<Real>& avg_jump2   = m_avgJump[a_dit].get_ivfab(phase2);
  
  const EBCellFAB& phi1              = a_phi.getPhase(phase1);
  const EBCellFAB& phi2              = a_phi.getPhase(phase2);

  const EBISBox& ebisbox1            = phi1.getEBISBox();
  const EBISBox& ebisbox2            = phi2.getEBISBox();
  
  BaseIVFAB<Real>& inhomo1 = m_inhomo[a_dit].get_ivfab(phase1);
  BaseIVFAB<Real>& inhomo2 = m_inhomo[a_dit].get_ivfab(phase2);

  BaseIVFAB<Real>& homog1  = m_homog[a_dit].get_ivfab(phase1);
  BaseIVFAB<Real>& homog2  = m_homog[a_dit].get_ivfab(phase2);

  // homog1.setVal(0.0); // This is the contribution of the surface charge for phase1
  // homog2.setVal(0.0); // This is the contribution of the surface charge for phase2

  // inhomo1.setVal(0.0); // This is the contribution from the other side for phase1
  // inhomo2.setVal(0.0); // This is the contribution from the other side for phase2
  
  // Compute the average jump. 
  for (IVSIterator ivsit(m_ivs[a_dit]); ivsit.ok(); ++ivsit){
    const IntVect iv = ivsit();

    const VolIndex vof0 = VolIndex(iv, 0);

    // Apply the stencils
    Real apply_sten1 = 0.0; // This is the applied stencil on phase 1. It will be an EB flux contribution to phase 2.
    Real apply_sten2 = 0.0; // This is the applied stencil on phase 2. It will be an EB flux contribution to phase 1.

    const VoFStencil& avgStencil1 = avg_sten1(vof0, 0);
    const VoFStencil& avgStencil2 = avg_sten2(vof0, 0);

    // Apply stencil on the gas side
    for (int i = 0; i < avgStencil1.size(); i++){
      const VolIndex ivof    = avgStencil1.vof(i);
      const Real     iweight = avgStencil1.weight(i);
      apply_sten1 += iweight*phi1(ivof,comp);
#if DEBUG
      const IntVect iiv = ivof.gridIndex();
      if(iweight != iweight)      MayDay::Abort("got weight2 NaN");
      if(ebisbox1.isCovered(iiv)) MayDay::Abort("reaching into covered cell for phase1");
#endif
    }

    // Apply stencil on the other side
    for (int i = 0; i < avgStencil2.size(); i++){
      const VolIndex ivof    = avgStencil2.vof(i);
      const Real     iweight = avgStencil2.weight(i);
      apply_sten2 += iweight*phi2(ivof,comp);
#if DEBUG
      const IntVect iiv = ivof.gridIndex();
      if(iweight != iweight)      MayDay::Abort("got weight2 NaN");
      if(ebisbox2.isCovered(iiv)) MayDay::Abort("reaching into covered cell for phase2");
#endif
    }

    // Compute the factor 1/(bp*wb + bq*wq)
    const Real factor = 1.0/(avg_bco1(vof0,0)*avg_w1(vof0,0) + avg_bco2(vof0,0)*avg_w2(vof0,0));

#if DEBUG 
    if(factor      != factor)      MayDay::Abort("got factor NaN");
    if(apply_sten1 != apply_sten1) MayDay::Abort("got sten1 NaN");
    if(apply_sten2 != apply_sten2) MayDay::Abort("got sten2 NaN");
#endif

    // Now make the "inhomogeneous" and "homogeneous" contributions
    const Vector<VolIndex> vofs1 = ebisbox1.getVoFs(iv);
    const Vector<VolIndex> vofs2 = ebisbox2.getVoFs(iv);

    // Phase 1 -- should the contributions be weighted somehow...?
    for (int v = 0; v < vofs1.size(); v++){
      const VolIndex vof = vofs1[v];

      homog1(vof,comp)  =  avg_jump1(vof0, 0);
      inhomo1(vof,comp) = -avg_bco2(vof0,0)*apply_sten2;

      homog1(vof,comp)  *= factor;
      inhomo1(vof,comp) *= factor;

#if DEBUG
      if(homog1(vof,comp)  != homog1(vof,comp))  MayDay::Abort("got homo1 NaN");
      if(inhomo1(vof,comp) != inhomo1(vof,comp)) MayDay::Abort("got inhomo1 NaN");
#endif
    }

    // Phase 2 -- should the contributions be weighted somehow...?
    for (int v = 0; v < vofs2.size(); v++){
      const VolIndex vof = vofs2[v];

      homog2(vof,comp)  =  avg_jump2(vof0, 0);
      inhomo2(vof,comp) = -avg_bco1(vof0,0)*apply_sten1;

      homog2(vof,comp)  *= factor;
      inhomo2(vof,comp) *= factor;

#if DEBUG
      if(homog2(vof,comp)  != homog2(vof,comp))  MayDay::Abort("got homo2 NaN");
      if(inhomo2(vof,comp) != inhomo2(vof,comp)) MayDay::Abort("got inhomo2 NaN");
#endif
    }
  }
}



#endif
