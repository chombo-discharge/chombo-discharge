#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NonAggregatedEBStencil.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "NamespaceHeader.H"
/**************/
/**************/
/**************/
NonAggregatedEBStencil::NonAggregatedEBStencil(const Vector<VolIndex>& a_srcVofs,
                                               const BaseIVFAB<VoFStencil>& a_vofStencil,
                                               const Box& a_box,
                                               const EBISBox& a_ebisBox,
                                               const IntVect& a_ghostVectPhi,
                                               const IntVect& a_ghostVectLph,
                                               int a_varDest,
                                               bool a_doRelaxOpt,
                                               int a_nComp, 
                                               IntVectSet a_setIrreg,
                                               bool       a_useInputSet)
{
  CH_TIME("NonAggregatedEBStencil::constructor1");
  m_srcVofs      = a_srcVofs;
  m_ebisBox      = a_ebisBox;
  m_destVar      = a_varDest;
  m_doRelaxOpt   = a_doRelaxOpt;
  m_nComp        = a_nComp;
  m_setIrreg     = a_setIrreg;
  m_useInputSets = a_useInputSet;
  m_cacheLph.resize(m_srcVofs.size());
  m_cachePhi.resize(m_srcVofs.size());
  m_vofStencil.resize(m_srcVofs.size());
  for (int isrc = 0; isrc < a_srcVofs.size(); isrc++)
    {
      m_vofStencil[isrc] = a_vofStencil(a_srcVofs[isrc], 0);
    }
}
/**************/
/**************/
NonAggregatedEBStencil::NonAggregatedEBStencil(const Vector<VolIndex>& a_srcVofs,
                                               const Vector<VoFStencil>& a_vofStencil,
                                               const Box& a_boxLph,
                                               const Box& a_boxPhi,
                                               const EBISBox& a_ebisBoxLph,
                                               const EBISBox& a_ebisBoxPhi,
                                               const IntVect& a_ghostVectLph,
                                               const IntVect& a_ghostVectPhi,
                                               int a_varDest,
                                               int a_nComp,
                                               IntVectSet a_setIrreg,
                                               bool       a_useInputSet)
{
  CH_TIME("NonAggregatedEBStencil::constructor2");
  m_srcVofs      = a_srcVofs;
  m_vofStencil   = a_vofStencil;
  m_ebisBox      = a_ebisBoxPhi;
  m_destVar      = a_varDest;
  m_nComp        = a_nComp;
  m_setIrreg     = a_setIrreg;
  m_useInputSets = a_useInputSet;
  m_cacheLph.resize(a_srcVofs.size());
  m_cachePhi.resize(a_srcVofs.size());
}
/**************/
/**************/

void NonAggregatedEBStencil::apply(EBCellFAB& a_lofphi, const EBCellFAB& a_phi, bool a_incrementOnly, int  a_ivar) const
{

  CH_TIME("NonAggregatedEBStencil::apply");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      Real& lphi = a_lofphi(vof, m_destVar);

      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      Real stenval = applyVoFStencil(m_vofStencil[isrc], a_phi, 0);
      lphi += stenval;
    }
}

void NonAggregatedEBStencil::apply(EBCellFAB&             a_lofphi,
                                   const EBCellFAB&       a_phi,
                                   const BaseIVFAB<Real>& a_alphaWeight,
                                   Real                   a_alpha,
                                   Real                   a_beta,
                                   bool                   a_incrementOnly) const

{
  CH_TIME("NonAggregatedEBStencil::apply1");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      Real& lphi = a_lofphi(vof, m_destVar);

      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      Real stenval   = applyVoFStencil(m_vofStencil[isrc], a_phi, 0);
      Real totalval  = a_alpha*a_alphaWeight(vof, 0)*a_phi(vof,0) + a_beta*stenval;
      lphi += totalval;
    }
}

//For EB x domain where m_alpha, m_beta have changed since defineStencils
//we need both alphaWeight and betaWeight to calculate the relaxation parameter:
//lambdaDiagWeight = 1/(alpha*alphaWeight+beta*betaWeight)
void NonAggregatedEBStencil::apply(EBCellFAB&             a_lofphi,
                                   const EBCellFAB&       a_phi,
                                   const Real             a_lambdaFactor,
                                   const Real             a_alpha,
                                   const BaseIVFAB<Real>& a_alphaWeight,
                                   const Real             a_beta,
                                   const BaseIVFAB<Real>& a_betaWeight,
                                   Real                   a_one,
                                   bool                   a_incrementOnly) const

{
  CH_TIME("NonAggregatedEBStencil::apply2");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      Real& lphi = a_lofphi(vof, m_destVar);

      if (!a_incrementOnly)
        {
          lphi =  0.;
        }
      //for some reason this version does not include lphi in EBStencil ---dtg
      //Real stenval   = applyVoFStencil(m_vofStencil[isrc], a_phi, 0);
      //      lphi += stenval
      Real phival    = a_phi(vof, 0);

      Real alphaWeight = a_alphaWeight(vof, 0);
      Real betaWeight =  a_betaWeight(vof, 0);
      Real product = a_alpha*alphaWeight+a_beta*betaWeight;
      Real lambdaWeight = 0.;
      if (Abs(product) > 1.e-15) lambdaWeight = 1./product;

      lphi += a_lambdaFactor*lambdaWeight*phival;
    }
}

void NonAggregatedEBStencil::applyInhomDomBC(EBCellFAB&             a_lofphi,
                                             const EBCellFAB&       a_phi,
                                             const Real             a_factor) const

{
  CH_TIME("NonAggregatedEBStencil::applyInhomDomBC");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      Real& lphi       = a_lofphi(vof, m_destVar);
      const Real& sour = a_phi(vof,0);
      lphi += a_factor*sour;
    }
}

void NonAggregatedEBStencil::relax(EBCellFAB&             a_phi,
                                   const EBCellFAB&       a_rhs,
                                   const BaseIVFAB<Real>& a_alphaWeight,
                                   const BaseIVFAB<Real>& a_betaWeight,
                                   Real                   a_alpha,
                                   Real                   a_beta,
                                   Real                   a_safety) const
{
  CH_TIME("NonAggregatedEBStencil::relax");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      Real& phi = a_phi(vof,0);
      Real rhs = a_rhs(vof, 0);
      Real stenval   = applyVoFStencil(m_vofStencil[isrc], a_phi, 0);
      Real alphaWeight = a_alphaWeight(vof, 0);
      Real betaWeight = a_betaWeight(vof, 0);
      Real denom = a_alpha*alphaWeight + a_beta*betaWeight;
      Real lambda = 0;
      if (Abs(denom) > 1.0e-12)
        lambda = a_safety/(denom);

      Real lphi = a_alpha*alphaWeight*phi + a_beta*stenval;

      phi = phi + lambda * (rhs - lphi);
    }
}

/**************/
void NonAggregatedEBStencil::relaxClone(EBCellFAB&             a_phi,
                                        const EBCellFAB&       a_phiOld,
                                        const EBCellFAB&       a_rhs,
                                        const BaseIVFAB<Real>& a_alphaWeight,
                                        const BaseIVFAB<Real>& a_betaWeight,
                                        Real a_alpha, Real a_beta, Real a_safety) const
{
  CH_TIME("NonAggregatedEBStencil::relaxClone");

  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      Real& phi = a_phi(vof,0);
      Real phiOld = a_phi(vof, 0);
      Real rhs = a_rhs(vof, 0);
      Real stenval   = applyVoFStencil(m_vofStencil[isrc], a_phi, 0);
      Real alphaWeight = a_alphaWeight(vof, 0);
      Real betaWeight = a_betaWeight(vof, 0);
      Real denom = a_alpha*alphaWeight + a_beta*betaWeight;
      Real lambda = 0;
      if (Abs(denom) > 1.0e-12)
        lambda = a_safety/(denom);

      Real lphi = a_alpha*alphaWeight*phi + a_beta*stenval;

      phi = phiOld + lambda * (rhs - lphi);
    }
}

/**************/
void
NonAggregatedEBStencil::cachePhi(const EBCellFAB& a_phi, int a_ivar) const
{
  CH_TIME("NonAggregatedEBStencil::cachePhi");

  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      m_cachePhi[isrc] = a_phi(vof, a_ivar);
    }
}

void
NonAggregatedEBStencil::cache(const EBCellFAB& a_lph, int a_ivar) const
{
  CH_TIME("NonAggregatedEBStencil::cacheLph");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      m_cacheLph[isrc] = a_lph(vof, a_ivar);
    }
}
/**************/
/**************/
void
NonAggregatedEBStencil::uncachePhi(EBCellFAB& a_phi, int a_ivar) const
{
  CH_TIME("NonAggregatedEBStencil::uncachePhi");

  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      a_phi(vof, a_ivar) = m_cachePhi[isrc];
    }
}
/**************/
void
NonAggregatedEBStencil::uncache(EBCellFAB& a_lph, int a_ivar) const
{
  CH_TIME("NonAggregatedEBStencil::uncacheLph");
  for (int isrc = 0; isrc < m_srcVofs.size(); isrc++)
    {
      const VolIndex& vof = m_srcVofs[isrc];
      a_lph(vof, a_ivar)  =   m_cacheLph[isrc];
    }
}

/**************/
/**************/
NonAggregatedEBStencil::~NonAggregatedEBStencil()
{
}
/**************/

#include "NamespaceFooter.H"
