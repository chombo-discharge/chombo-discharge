#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "dirichletconductivityebbc.H"

#include "BoxIterator.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"

namespace ChomboDischarge {

  void
  dirichletconductivityebbc::
  define(const LayoutData<IntVectSet>& a_cfivs,
	 const Real&                   a_factor)
  {
    if (!m_coefSet)
      {
	MayDay::Error("DirCondEBBC: need to call setCoef BEFORE calling define.");
      }
    m_bc.define(a_cfivs, a_factor);
    LayoutData<BaseIVFAB<VoFStencil> >& poissSten = *(m_bc.getFluxStencil(0));
    BoxLayout dbl = poissSten.boxLayout();
    m_fluxStencil.define(dbl);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
	const IntVectSet& ivs = poissSten[dit()].getIVS();
	const EBGraph&    ebg = poissSten[dit()].getEBGraph();
	m_fluxStencil[dit()].define(ivs, ebg, 1);
	for (VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
	  {
	    m_fluxStencil[dit()](vofit(), 0) = poissSten[dit()](vofit(), 0);
	    Real factor = (*m_bcoe)[dit()](vofit(), 0);
	    factor *= m_beta;
	    m_fluxStencil[dit()](vofit(), 0) *= factor;
	  }
      }
  }
  /*****************/
  void
  dirichletconductivityebbc::
  getEBFlux(Real&                         a_flux,
	    const VolIndex&               a_vof,
	    const LevelData<EBCellFAB>&   a_phi,
	    const LayoutData<IntVectSet>& a_cfivs,
	    const DataIndex&              a_dit,
	    const RealVect&               a_probLo,
	    const RealVect&               a_dx,
	    const bool&                   a_useHomogeneous,
	    const Real&                   a_time,
	    const pair<int,Real>*         a_cacheHint )
  {
#if 1 // Robert, change Dec. 20, 2017. This can be removed.
    MayDay::Abort("how did I get called?");
#endif
    m_bc.getEBFlux(a_flux,
		   a_vof,
		   a_phi,
		   a_cfivs,
		   a_dit,
		   a_probLo,
		   a_dx,
		   a_useHomogeneous,
		   a_time,
		   a_cacheHint );

    Real bcoef = (*m_bcoe)[a_dit](a_vof,0);
    a_flux *= bcoef;
  }
  /*****************/
  dirichletconductivityebbc::
  ~dirichletconductivityebbc()
  {
  }
  /*****************/
  void
  dirichletconductivityebbc::
  setValue(Real a_value)
  {
    m_dataBased = false;
    m_bc.setValue(a_value);
  }
  /*****************/
  void
  dirichletconductivityebbc::
  setFunction(RefCountedPtr<BaseBCValue> a_flux)
  {
    m_dataBased = false;
    m_bc.setFunction(a_flux);
  }
  /*****************/
  void
  dirichletconductivityebbc::
  applyEBFlux(EBCellFAB&                    a_lphi,
	      const EBCellFAB&              a_phi,
	      VoFIterator&                  a_vofit,
	      const LayoutData<IntVectSet>& a_cfivs,
	      const DataIndex&              a_dit,
	      const RealVect&               a_probLo,
	      const RealVect&               a_dx,
	      const Real&                   a_factor,
	      const bool&                   a_useHomogeneous,
	      const Real&                   a_time)
  {
    CH_TIME("dirichletconductivityebbc::applyEBFlux");
    CH_assert(a_lphi.nComp() == 1 );
    CH_assert(a_phi.nComp()  == 1);

    if(!a_useHomogeneous){
      const BaseIVFAB<Real>& poissWeight = (m_bc.getFluxWeight())[a_dit];
      Real value = 0.0;

      const EBISBox&   ebisBox = a_phi.getEBISBox();
      for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
	{
	  const VolIndex& vof = a_vofit();
	  if (m_dataBased)
	    {
	      //          if ((*m_data)[a_dit].getIVS().contains(vof.gridIndex()))
	      //             {
	      //               value = (*m_data)[a_dit](vof, 0);
	      //             }
	      //          else
	      //            {
	      //              value = 0.;
	      //            }
	      value = (*m_data)[a_dit](vof, 0);

	    }
	  else if (m_bc.m_isFunction)
	    {
	      const RealVect& centroid = ebisBox.bndryCentroid(vof);
	      const RealVect&   normal = ebisBox.normal(vof);

	      value = m_bc.m_func->value(vof,centroid,normal,a_dx,a_probLo,a_dit,a_time,0);
	    }
	  else
	    {
	      if (m_bc.m_onlyHomogeneous)
		{
		  MayDay::Error("dirichletconductivityebbc::getFaceFlux called with undefined inhomogeneous BC");
		}

	      value = m_bc.m_value;
	    }
	  Real poissWeightPt = poissWeight(vof, 0);
	  const Real& areaFrac = ebisBox.bndryArea(vof);
	  const Real& bcoef    = (*m_bcoe)[a_dit](vof,0);
	  Real flux = poissWeightPt*value*areaFrac;
	  Real compFactor = a_factor*bcoef*m_beta;
	  a_lphi(vof,0) += flux * compFactor;
	}
    }
  }
  /*****************/
  dirichletconductivityebbcfactory::
  dirichletconductivityebbcfactory()
  {
    m_value = 12345.6789;
    m_flux = RefCountedPtr<BaseBCValue>();
    m_onlyHomogeneous = true;
    m_isFunction = false;
    m_dataBased = false;
  }
  /*****************/
  dirichletconductivityebbcfactory::
  ~dirichletconductivityebbcfactory()
  {
  }
  /*****************/
  void
  dirichletconductivityebbcfactory::
  setValue(Real a_value)
  {
    m_value = a_value;
    m_flux = RefCountedPtr<BaseBCValue>();

    m_onlyHomogeneous = false;
    m_isFunction = false;
  }
  /*****************/
  void
  dirichletconductivityebbcfactory::
  setFunction(RefCountedPtr<BaseBCValue> a_flux)
  {
    m_value = 12345.6789;
    m_flux = a_flux;

    m_onlyHomogeneous = false;
    m_isFunction = true;
  }
  /*****************/
  dirichletconductivityebbc*
  dirichletconductivityebbcfactory::
  create(const ProblemDomain& a_domain,
	 const EBISLayout&    a_layout,
	 const RealVect&      a_dx,
	 const IntVect*       a_ghostCellsPhi /*=0*/,
	 const IntVect*       a_ghostCellsRhs /*=0*/)
  {
    dirichletconductivityebbc* fresh = new dirichletconductivityebbc(a_domain,a_layout,a_dx, a_ghostCellsPhi, a_ghostCellsRhs);
    fresh->setOrder(m_order);

    if (!m_onlyHomogeneous)
      {
	if (m_dataBased)
	  {
	    fresh->setData(m_data);
	  }
	else if (!m_isFunction)
	  {
	    fresh->setValue(m_value);
	  }
	else
	  {
	    fresh->setFunction(m_flux);
	  }
      }

    return fresh;
  }
}
#include "NamespaceFooter.H"
