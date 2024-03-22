/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFLevelGrid.cpp
  @brief  Implementation of CD_MFLevelGrid.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_MFLevelGrid.H>
#include <CD_NamespaceHeader.H>

MFLevelGrid::MFLevelGrid()
{}

void
MFLevelGrid::setMaxRefinementRatio(const int a_refRat)
{
  for (int i = 0; i < m_eblg.size(); i++) {
    m_eblg[i].setMaxRefinementRatio(a_refRat);
  }
}

// Olde code below here

MFLevelGrid::MFLevelGrid(const DisjointBoxLayout&                   a_dbl,
                         const ProblemDomain&                       a_domain,
                         const int                                  a_ebghost,
                         const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace)
{
  m_multifluidIndexSpace = a_multiFluidIndexSpace;
  m_eblg.resize(0);
  for (int i = 0; i < a_multiFluidIndexSpace->numPhases(); i++) {
    m_eblg.push_back(EBLevelGrid(a_dbl, a_domain, a_ebghost, a_multiFluidIndexSpace->getEBIndexSpace(i)));
  }
}

MFLevelGrid::MFLevelGrid(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace,
                         const Vector<EBLevelGrid>&                 a_eblg)
{
  this->define(a_multiFluidIndexSpace, a_eblg);
}

MFLevelGrid::~MFLevelGrid()
{}

int
MFLevelGrid::numPhases() const
{
  return m_eblg.size();
}

void
MFLevelGrid::define(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace,
                    const Vector<EBLevelGrid>&                 a_eblg)
{
  m_multifluidIndexSpace = a_multiFluidIndexSpace;
  m_eblg                 = a_eblg;
}

const RefCountedPtr<MultiFluidIndexSpace>&
MFLevelGrid::getMfIndexSpace() const
{
  return m_multifluidIndexSpace;
}

ProblemDomain
MFLevelGrid::getDomain() const
{
  return m_eblg[0].getDomain();
}

DisjointBoxLayout
MFLevelGrid::getGrids() const
{
  return m_eblg[0].getDBL();
}

EBLevelGrid&
MFLevelGrid::getEBLevelGrid(int a_phase)
{
  CH_assert(a_phase < m_eblg.size());
  return m_eblg[a_phase];
}

const EBLevelGrid&
MFLevelGrid::getEBLevelGrid(int a_phase) const
{
  CH_assert(a_phase < m_eblg.size());
  return m_eblg[a_phase];
}

IntVectSet
MFLevelGrid::interfaceRegion(const Box& a_box, const DataIndex& a_dit, const int a_phase1, const int a_phase2) const
{

  IntVectSet ret;

  if (m_multifluidIndexSpace->numPhases() == 2) {
    const EBLevelGrid& eblg1 = this->getEBLevelGrid(a_phase1);
    const EBLevelGrid& eblg2 = this->getEBLevelGrid(a_phase2);

    const EBISBox& ebisbox1 = eblg1.getEBISL()[a_dit];
    const EBISBox& ebisbox2 = eblg2.getEBISL()[a_dit];

    ret = ebisbox1.getIrregIVS(a_box) & ebisbox2.getIrregIVS(a_box);
  }

  return ret;
}

#if 0
bool MFLevelGrid::interfacePair(IntVect&             a_iv,
				 const IntVect&       a_iv_in,
				 const Box&           a_box,
				 const DataIndex&     a_dit,
				 const int            a_phase1,
				 const int            a_phase2) const {
  bool found_iv = false;

  if(m_multifluidIndexSpace->numPhases() == 2){
    const EBLevelGrid& eblg1 = this->getEBLevelGrid(a_phase1);
    const EBLevelGrid& eblg2 = this->getEBLevelGrid(a_phase2);
      
    const EBISBox& ebisbox1  = eblg1.getEBISL()[a_dit];
    const EBISBox& ebisbox2  = eblg2.getEBISL()[a_dit];

    const IntVectSet& irreg2 = ebisbox2.getIrregIVS(a_box);

    // Pairing IV must be the current iv, or one of its four (six) neighbors
    if(irreg2.contains(a_iv_in)){
      a_iv     = a_iv_in;
      found_iv = true;
    }
    else{
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  const IntVect iv = a_iv_in + BASISV(dir)*sign(sit());
	  if(irreg2.contains(iv)){
	    if(found_iv){
	      MayDay::Abort("MFLevelGrid::interfacePair - iv has multiple neighboring irregular cells. Aborting.");
	    }
	    else{
	      a_iv     = iv;
	      found_iv = true;
	    }
	  }
	}
      }
    }
  }

  return found_iv;
}
#endif

#include <CD_NamespaceFooter.H>
