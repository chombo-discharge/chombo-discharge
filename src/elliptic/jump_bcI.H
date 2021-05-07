/*!
  @file   jump_bcI.H
  @brief  Implementation of jump_bc.cpp
  @author Robert Marskar
  @date   Dec. 2017
*/

#ifndef _JUMP_BCI_
#define _JUMP_BCI_

#include "jump_bc.H"

namespace ChomboDischarge {

#define DEBUG_JUMP 0

  inline
  void jump_bc::compute_avg_jump(const BaseIVFAB<Real>& a_jump, const MFCellFAB& a_phi, const DataIndex& a_dit){

    const int comp   = 0;
    const int phase1 = 0;
    const int phase2 = 1;
  
    const BaseIVFAB<Real>& bco1 = m_bco[a_dit].get_ivfab(phase1);
    const BaseIVFAB<Real>& bco2 = m_bco[a_dit].get_ivfab(phase2);

    const EBISBox& ebisbox1 = a_phi.getPhase(phase1).getEBISBox();
    const EBISBox& ebisbox2 = a_phi.getPhase(phase2).getEBISBox();

    BaseIVFAB<Real>& avg1 = m_avgJump[a_dit].get_ivfab(phase1);
    BaseIVFAB<Real>& avg2 = m_avgJump[a_dit].get_ivfab(phase2);

    avg1.setVal(0.0);
    avg2.setVal(0.0);

    for (IVSIterator ivsit(m_ivs[a_dit]); ivsit.ok(); ++ivsit){
      const IntVect iv = ivsit();

      const VolIndex vof0(iv, comp);
    
      const Vector<VolIndex> vofs1 = ebisbox1.getVoFs(iv);
      const Vector<VolIndex> vofs2 = ebisbox2.getVoFs(iv);

      // Compute average jump condition. a_jump should be defined on phase 1. 
      Real totalJump = 0.0;
      for (int ivof = 0; ivof < vofs1.size(); ivof++){
	const VolIndex& vof = vofs1[ivof];

	totalJump += a_jump(vof, comp);
      }
      totalJump /= vofs1.size();


      // Set the jump. It should not matter if the cell is multivalued because
      // the other VoFs in the cells should never be referenced. 
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

    jump_bc::compute_avg_jump(a_jump, a_phi, a_dit);
  
    const BaseIVFAB<Real>& avg_bco1    = m_avgBco[a_dit].get_ivfab(phase1);
    const BaseIVFAB<Real>& avg_bco2    = m_avgBco[a_dit].get_ivfab(phase2);
  
    const BaseIVFAB<Real>& avg_w1      = m_avgWeights[a_dit].get_ivfab(phase1);
    const BaseIVFAB<Real>& avg_w2      = m_avgWeights[a_dit].get_ivfab(phase2);
  
    const BaseIVFAB<VoFStencil>& avg_sten1 = m_avgStencils[a_dit].get_ivfab(phase1);
    const BaseIVFAB<VoFStencil>& avg_sten2 = m_avgStencils[a_dit].get_ivfab(phase2);

    const BaseIVFAB<Real>& avg_jump1   = m_avgJump[a_dit].get_ivfab(phase1);
    const BaseIVFAB<Real>& avg_jump2   = m_avgJump[a_dit].get_ivfab(phase2);

    const BaseIVFAB<Real>& avg_factor1   = m_avgFactor[a_dit].get_ivfab(phase1);
    const BaseIVFAB<Real>& avg_factor2   = m_avgFactor[a_dit].get_ivfab(phase2);
  
    const EBCellFAB& phi1              = a_phi.getPhase(phase1);
    const EBCellFAB& phi2              = a_phi.getPhase(phase2);

    const EBISBox& ebisbox1            = phi1.getEBISBox();
    const EBISBox& ebisbox2            = phi2.getEBISBox();
  
    BaseIVFAB<Real>& inhomo1 = m_inhomo[a_dit].get_ivfab(phase1);
    BaseIVFAB<Real>& inhomo2 = m_inhomo[a_dit].get_ivfab(phase2);

    BaseIVFAB<Real>& homog1  = m_homog[a_dit].get_ivfab(phase1);
    BaseIVFAB<Real>& homog2  = m_homog[a_dit].get_ivfab(phase2);

    homog1.setVal(0.0); // This is the contribution of the surface charge for phase1
    homog2.setVal(0.0); // This is the contribution of the surface charge for phase2

    inhomo1.setVal(0.0); // This is the contribution from the other side for phase1
    inhomo2.setVal(0.0); // This is the contribution from the other side for phase2
  
    // Compute the average jump. 
    for (IVSIterator ivsit(m_ivs[a_dit]); ivsit.ok(); ++ivsit){
      const IntVect iv = ivsit();

      const VolIndex vof0 = VolIndex(iv, 0);

      const Real& jump1 = avg_jump1(vof0, 0);
      const Real& jump2 = avg_jump2(vof0, 0);

      const Real& bco1  = avg_bco1(vof0, 0);
      const Real& bco2  = avg_bco2(vof0, 0);

      const Real& w1    = avg_w1(vof0, 0);
      const Real& w2    = avg_w2(vof0, 0);

      const VoFStencil& avgStencil1 = avg_sten1(vof0, 0);
      const VoFStencil& avgStencil2 = avg_sten2(vof0, 0);

      //    const Real factor = 1.0/(bco1*w1 + bco2*w2);
      const Real& factor1 = avg_factor1(vof0, 0);
      const Real& factor2 = avg_factor2(vof0, 0);

      const Vector<VolIndex> vofs1 = ebisbox1.getVoFs(iv);
      const Vector<VolIndex> vofs2 = ebisbox2.getVoFs(iv);

#if DEBUG_JUMP
      if (ebisbox1.isCovered(iv)) MayDay::Abort("jump_bc::match_bc -- cell is covered on phase 1");
      if (ebisbox2.isCovered(iv)) MayDay::Abort("jump_bc::match_bc -- cell is covered on phase 2");
#endif

      // Apply the stencils on each side. apply_sten1 is an EB flux contribution to phase 2 and vice versa. 
      Real apply_sten1 = 0.0;
      Real apply_sten2 = 0.0;

      for (int i = 0; i < avgStencil1.size(); i++){
	const VolIndex& ivof    = avgStencil1.vof(i);
	const Real&     iweight = avgStencil1.weight(i);
	apply_sten1 += iweight*phi1(ivof,comp);
      }

      for (int i = 0; i < avgStencil2.size(); i++){
	const VolIndex& ivof    = avgStencil2.vof(i);
	const Real&     iweight = avgStencil2.weight(i);
	apply_sten2 += iweight*phi2(ivof,comp);
      }

#if DEBUG_JUMP 
      if(std::isnan(apply_sten1)) MayDay::Abort("got sten1 NaN");
      if(std::isnan(apply_sten2)) MayDay::Abort("got sten2 NaN");
#endif

      // Now make the "inhomogeneous" and "homogeneous" contributions on each side. This gets area weighted in mfdirichletconductivityebbc. 
      for (int v = 0; v < vofs1.size(); v++){
	const VolIndex vof = vofs1[v];

	homog1(vof,comp)  =  jump1;
	inhomo1(vof,comp) = -bco2*apply_sten2;

	homog1(vof,comp)  *= factor1;
	inhomo1(vof,comp) *= factor1;

#if DEBUG_JUMP
	if(std::isnan(homog1(vof,comp)))  MayDay::Abort("got homo1 NaN");
	if(std::isnan(inhomo1(vof,comp))) MayDay::Abort("got inhomo1 NaN");
#endif
      }

      for (int v = 0; v < vofs2.size(); v++){
	const VolIndex vof = vofs2[v];

	homog2(vof,comp)  =  jump2;
	inhomo2(vof,comp) = -bco1*apply_sten1; 

	homog2(vof,comp)  *= factor2;
	inhomo2(vof,comp) *= factor2;

#if DEBUG_JUMP
	if(std::isnan(homog2(vof,comp)))  MayDay::Abort("got homo2 NaN");
	if(std::isnan(inhomo2(vof,comp))) MayDay::Abort("got inhomo2 NaN");
#endif
      }
    }
  }
}


#endif
