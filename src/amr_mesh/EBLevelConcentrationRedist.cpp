/*!
  @file   EBLevelConcentrationRedist.cpp
  @brief  Declaration file for EBLevelConcentrationRedist
  @author Robert Marskar
*/

#include "EBLevelConcentrationRedist.H"
#include "BaseIVFactory.H"

EBLevelConcentrationRedist::EBLevelConcentrationRedist(){
  m_isDefined = false;
}


EBLevelConcentrationRedist::~EBLevelConcentrationRedist(){

}

EBLevelConcentrationRedist::EBLevelConcentrationRedist(const DisjointBoxLayout& a_dbl,
						       const EBISLayout&        a_ebisl,
						       const ProblemDomain&     a_domain,
						       const int                a_ncomp,
						       const int                a_redistRad){
  CH_TIME("EBLevelConcentrationRedist::EBLevelConcentrationRedist");
  define(a_dbl, a_ebisl, a_domain, a_ncomp, a_redistRad);
}


void EBLevelConcentrationRedist::define(const DisjointBoxLayout& a_dbl,
					const EBISLayout&        a_ebisl,
					const ProblemDomain&     a_domain,
					const int                a_ncomp,
					const int                a_redistRad){
  CH_TIME("EBLevelConcentrationRedist::define");
  m_isDefined = true;
  m_grids     = a_dbl;
  m_ebisl     = a_ebisl;
  m_domain    = a_domain;
  m_redistRad = a_redistRad;
  m_ncomp    = a_ncomp;

  m_stencil.define(m_grids, m_ebisl, m_domain, m_redistRad, false);

  //define the sets for iterating over
  //to be the irregular cells over the
  //grid grown by the redistribution radius
  m_sets.define(m_grids);
  int redistRad = m_stencil.getRedistRadius();

  DataIterator dit = m_grids.dataIterator(); 
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      Box thisBox = m_grids.get(dit[mybox]);
      thisBox.grow(redistRad);
      thisBox &= m_domain;
      m_sets[dit[mybox]] = m_ebisl[dit[mybox]].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_ebisl, m_sets);
  IntVect ivghost = redistRad*IntVect::Unit;
  m_buffer.define(m_grids, m_ncomp, ivghost, factory);
  setToZero();
  
}

void EBLevelConcentrationRedist::resetWeights(const LevelData<EBCellFAB>& a_modifier, const int a_ivar){
  CH_TIME("EBLevelConcentrationRedist::resetWeights");
}

void EBLevelConcentrationRedist::increment(const BaseIVFAB<Real>& a_massDiff,
					   const DataIndex&       a_datInd,
					   const Interval&        a_variables){
  CH_TIME("EBLevelConcentrationRedist::increment");

  CH_assert(m_isDefined);
  BaseIVFAB<Real>& bufFAB = m_buffer[a_datInd];
  const IntVectSet& fabIVS = a_massDiff.getIVS();
  const IntVectSet& bufIVS = m_sets[a_datInd];

  IntVectSet ivs = m_ebisl[a_datInd].getIrregIVS(m_grids.get(a_datInd));;
  CH_assert(fabIVS.contains(ivs));
  CH_assert(bufIVS.contains(ivs));
  for (VoFIterator vofit(ivs, m_ebisl[a_datInd].getEBGraph());
       vofit.ok(); ++vofit)
    {
      for (int ivar = a_variables.begin();
	   ivar <= a_variables.end(); ivar++)
        {
          bufFAB(vofit(), ivar) += a_massDiff(vofit(), ivar);
        }
    }
}

void EBLevelConcentrationRedist::redistribute(LevelData<EBCellFAB>& a_solution, const Interval& a_variables){
  CH_TIME("EBLevelConcentrationRedist::redistribute(LDEBCellFAB, Interval");

  this->redistribute(a_solution, a_variables, a_variables);
}
		 
void EBLevelConcentrationRedist::redistribute(LevelData<EBCellFAB>& a_solution,
					      const Interval& a_srcVar,
					      const Interval& a_dstVar){
  CH_TIME("EBLevelConcentrationRedist::redistribute(LDEBCellFAB, Interval, Interval");
  CH_assert(a_srcVar.size() == a_dstVar.size());
  CH_assert(a_srcVar.begin() >= 0);
  CH_assert(a_dstVar.begin() >= 0);
  CH_assert(a_srcVar.end() <  m_ncomp);
  CH_assert(a_dstVar.end() <  a_solution.nComp());

  CH_assert(m_isDefined);

  //pout() << "redistribute 0" << endl;
  //exchange ghost cell information of the buffer.
  //this way the redistribution from ghost cells will
  //account for fine-fine interfaces.
  Interval wholeInterv(0, m_ncomp-1);
  m_buffer.exchange(wholeInterv);
  //pout() << "redistribute 1" << endl;
  //loop over grids.
  DataIterator dit = m_grids.dataIterator(); 
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      const BaseIVFAB<Real>& bufFAB = m_buffer[dit[mybox]];
      const BaseIVFAB<VoFStencil>& stenFAB = m_stencil[dit[mybox]];
      const Box& grid = m_grids.get(dit[mybox]);
      EBCellFAB& solFAB = a_solution[dit[mybox]];

      //pout() << "redistribute 2" << endl;
      //loop over the irregular vofs in the grown box.
      for (VoFIterator vofit(m_sets[dit[mybox]], m_ebisl[dit[mybox]].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& srcVoF = vofit();
          const VoFStencil& vofsten = stenFAB(srcVoF, 0);
          for (int isten = 0; isten < vofsten.size(); isten++)
            {
              const Real& weight = vofsten.weight(isten);
              const VolIndex& dstVoF = vofsten.vof(isten);
              //only add the mass to the solution if it is in the
              //solutions valid region --- because the buffer is over
              //the grown box, this prevents redistribution off
              //the end of the world.
              if (grid.contains(dstVoF.gridIndex()))
                {
                  for (int ivar = 0; ivar < a_srcVar.size(); ivar++)
                    {
                      int isrc = ivar + a_srcVar.begin();
                      int idst = ivar + a_dstVar.begin();

                      const Real& mass = bufFAB(srcVoF, isrc);
                      const Real& solu = solFAB(dstVoF, idst);
                      solFAB(dstVoF, idst) = mass*weight + solu;
                    }
                  //ch_flops()+=a_srcVar.size()*2; 
                } //end check if adding to valid soltuion
            } //end loop over vof stencil
        } //end loop over irregular vofs in grown box
    } //end loop over grids.

}

void EBLevelConcentrationRedist::setToZero() {
  CH_TIME("EBLevelConcentrationRedist::setToZero");
  CH_assert(m_isDefined);
  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    m_buffer[dit()].setVal(0.0);
  }
}
