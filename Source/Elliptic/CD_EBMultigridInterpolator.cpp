/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.cpp
  @brief  Implementation of CD_EBMultigridInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <EBLevelDataOps.H>
#include <NeighborIterator.H>
#include <InterpF_F.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

EBMultigridInterpolator::EBMultigridInterpolator(){
  m_isDefined = false;
}

EBMultigridInterpolator::EBMultigridInterpolator(const EBLevelGrid& a_eblgFine,
						 const EBLevelGrid& a_eblgCoar,
						 const int          a_nRef,
						 const int          a_nVar,
						 const int          a_ghostCF) {

  if(a_ghostCF != 1) MayDay::Abort("EBMultigridInterpolator::EBMultigridInterpolator - only one ghost cell supported (for now)!");
  if(a_ghostCF <= 0) MayDay::Abort("EBMultigridInterpolator::EBMultigridInterpolator - must interpolator at least one ghost cell!");
  
  m_ghostCF = a_ghostCF;

  // Build the CFIVS for EBQuadCFInterp. 
  const DisjointBoxLayout& dblFine = a_eblgFine.getDBL();
  const ProblemDomain& domainFine  = a_eblgFine.getDomain();
  
  LayoutData<IntVectSet> cfivs(dblFine);
  for (DataIterator dit(dblFine); dit.ok(); ++dit){
    Box grownBox = grow(dblFine[dit()], m_ghostCF);
    grownBox &= domainFine;

    cfivs[dit()] = IntVectSet(grownBox);
    
    NeighborIterator nit(dblFine);
    for (nit.begin(dit()); nit.ok(); ++ nit){
      cfivs[dit()] -= dblFine[nit()];
    }
  }

  EBQuadCFInterp::define(a_eblgFine.getDBL(),
			 a_eblgCoar.getDBL(),
			 a_eblgFine.getEBISL(),
			 a_eblgCoar.getEBISL(),
			 a_eblgCoar.getDomain(),
			 a_nRef,
			 a_nVar,
			 cfivs,
			 a_eblgFine.getEBIS(),
			 true);

  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;

  this->defineCFIVS();
  
  // Define a temp which is zero everywhere. Don't need ghost cells because the stencils should(!) only reach into valid cells AFAIK. 
  EBCellFactory cellFact(m_eblgCoar.getEBISL());
  m_zeroCoar.define(m_eblgCoar.getDBL(), a_nVar, IntVect::Zero, cellFact);
  EBLevelDataOps::setToZero(m_zeroCoar);
}

int EBMultigridInterpolator::getGhostCF() const{
  return m_ghostCF;
}

void EBMultigridInterpolator::defineCFIVS(){
  for (int dir = 0; dir < SpaceDim; dir++){
    m_loCFIVS[dir].define(m_eblgFine.getDBL());
    m_hiCFIVS[dir].define(m_eblgFine.getDBL());

    for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit){
      m_loCFIVS[dir][dit()].define(m_eblgFine.getDomain(), m_eblgFine.getDBL()[dit()], m_eblgFine.getDBL(), dir, Side::Lo);
      m_hiCFIVS[dir][dit()].define(m_eblgFine.getDomain(), m_eblgFine.getDBL()[dit()], m_eblgFine.getDBL(), dir, Side::Hi);
    }
  }
}

EBMultigridInterpolator::~EBMultigridInterpolator(){

}

void EBMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables){
  EBQuadCFInterp::interpolate(a_phiFine, a_phiCoar, a_variables);
}

void EBMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables){
#if 0 // This is very slow code
  EBQuadCFInterp::interpolate(a_phiFine, m_zeroCoar, a_variables);
#else // Much faster code
  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit){
    this->coarseFineInterpH(a_phiFine[dit()], a_variables, dit());
  }
#endif
}

void EBMultigridInterpolator::coarseFineInterpH(EBCellFAB& a_phi, const Interval a_variables, const DataIndex& a_dit){
  const Real m_dx     = 1.0;
  const Real m_dxCoar = m_dx*m_refRat;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){

      const CFIVS* cfivsPtr = nullptr;
      if(sit() == Side::Lo) {
	cfivsPtr = &m_loCFIVS[dir][a_dit];
      }
      else{
	cfivsPtr = &m_hiCFIVS[dir][a_dit];
      }

      const IntVectSet& ivs = cfivsPtr->getFineIVS();
      if (cfivsPtr->isPacked() ){
	const int ihiorlo = sign(sit());
	FORT_INTERPHOMO(CHF_FRA(a_phi.getSingleValuedFAB()),
			CHF_BOX(cfivsPtr->packedBox()),
			CHF_CONST_REAL(m_dx),
			CHF_CONST_REAL(m_dxCoar),
			CHF_CONST_INT(dir),
			CHF_CONST_INT(ihiorlo));
      }
      else {
	if(!ivs.isEmpty()){
	  
	  Real halfdxcoar = m_dxCoar/2.0;
	  Real halfdxfine = m_dx/2.0;
	  Real xg = halfdxcoar -   halfdxfine;
	  Real xc = halfdxcoar +   halfdxfine;
	  Real xf = halfdxcoar + 3*halfdxfine;
	  Real hf = m_dx;
	  Real denom = xf*xc*hf;

	  const EBISBox& ebisBox = m_eblgFine.getEBISL()[a_dit];
	  const EBGraph& ebgraph = ebisBox.getEBGraph();
	      
	  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	    const VolIndex& VoFGhost = vofit();

	    IntVect ivGhost  = VoFGhost.gridIndex();

	    Vector<VolIndex> farVoFs;
	    Vector<VolIndex> closeVoFs = ebisBox.getVoFs(VoFGhost, dir, flip(sit()), 1);
	    
	    bool hasClose = (closeVoFs.size() > 0);
	    bool hasFar = false;
	    Real phic = 0.0;
	    Real phif = 0.0;
	    
	    if (hasClose){
	      const int& numClose = closeVoFs.size();
	      for (int iVof=0;iVof<numClose;iVof++){
		const VolIndex& vofClose = closeVoFs[iVof];
		phic += a_phi(vofClose,0);
	      }
	      phic /= Real(numClose);

	      farVoFs = ebisBox.getVoFs(VoFGhost, dir, flip(sit()), 2);
	      hasFar   = (farVoFs.size()   > 0);
	      if (hasFar){
		const int& numFar = farVoFs.size();
		for (int iVof=0;iVof<numFar;iVof++){
		  const VolIndex& vofFar = farVoFs[iVof];
		  phif += a_phi(vofFar,0);
		}
		phif /= Real(numFar);
	      }
	    }

	    Real phiGhost;
	    if (hasClose && hasFar){
	      // quadratic interpolation  phi = ax^2 + bx + c
	      Real A = (phif*xc - phic*xf)/denom;
	      Real B = (phic*hf*xf - phif*xc*xc + phic*xf*xc)/denom;

	      phiGhost = A*xg*xg + B*xg;
	    }
	    else if (hasClose){
	      //linear interpolation
	      Real slope =  phic/xc;
	      phiGhost   =  slope*xg;
	    }
	    else{
	      phiGhost = 0.0; //nothing to interpolate from
	    }
	    a_phi(VoFGhost, 0) = phiGhost;
	  }
	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
