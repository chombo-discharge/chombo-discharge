/*!
  @file MFLevelGrid
  @brief Declaration of a wrapper for wrapping multifluid EBLevelGrids
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate header and implementation files
*/

#ifndef _MFLEVELGRID_
#define _MFLEVELGRID_

#include "mfis.H"

#include "EBLevelGrid.H"

/*!
  @brief Wrapper class for holding multifluid EBLevelGrids
*/
class MFLevelGrid {
public:

  /*!
    @brief Weak constructor, must subsequently call define
  */
  MFLevelGrid(){

  }

  /*!
    @brief Full constructor. Does full allocation
  */
  MFLevelGrid(const DisjointBoxLayout&          a_dbl,
	      const ProblemDomain&              a_domain,
	      const int                         a_ebghost,
	      const RefCountedPtr<mfis>&        a_mfis){
    m_mfis = a_mfis;
    m_eblg.resize(0);
    for (int i = 0; i < a_mfis->num_phases(); i++){
      m_eblg.push_back(EBLevelGrid(a_dbl, a_domain, a_ebghost, a_mfis->get_ebis(i)));
    }
  }

  /*!
    @brief Alternative full constructor - RefCounted storage so no allocation here
  */
  MFLevelGrid(const RefCountedPtr<mfis>& a_mfis,
	      const Vector<EBLevelGrid>& a_eblg){
    this->define(a_mfis, a_eblg);
  }

  /*!
    @brief Destructor
  */
  ~MFLevelGrid(){

  }

  /*!
    @brief Number of phases
  */
  virtual const int num_phases() const {
    return m_eblg.size();
  }

  /*!
    @brief Define function
  */
  virtual void define(const RefCountedPtr<mfis>& a_mfis,
		      const Vector<EBLevelGrid>& a_eblg){
    m_mfis = a_mfis;
    m_eblg = a_eblg;
  }

  /*!
    @brief Get the index space
  */
  virtual const RefCountedPtr<mfis>& get_mfis() const {
    return m_mfis;
  }

  /*!
    @brief Get the problem domain
  */
  virtual ProblemDomain get_domain() const {
    return m_eblg[0].getDomain();
  }

  /*!
    @brief Get the grids
  */
  virtual const DisjointBoxLayout get_grids() const {
    return m_eblg[0].getDBL();
  }

  /*!
    @brief Get phase
  */
  virtual EBLevelGrid& get_eblg(int a_phase){
    return m_eblg[a_phase];
  }

  /*!
    @brief Get phase
  */
  virtual const EBLevelGrid& get_eblg(int a_phase) const {
    return m_eblg[a_phase];
  }

  /*!
    @brief Get interface region between two phases.
  */
  virtual IntVectSet interface_region(const Box&       a_box,
				      const DataIndex& a_dit,
				      const int        a_phase1 = 0,
				      const int        a_phase2 = 1) const {

    IntVectSet ret;
    
    if(m_mfis->num_phases() == 2){
      const EBLevelGrid& eblg1 = this->get_eblg(a_phase1);
      const EBLevelGrid& eblg2 = this->get_eblg(a_phase2);
      
      const EBISBox& ebisbox1  = eblg1.getEBISL()[a_dit];
      const EBISBox& ebisbox2  = eblg2.getEBISL()[a_dit];

      ret = ebisbox1.getIrregIVS(a_box) & ebisbox2.getIrregIVS(a_box);
    }

    return ret;
  }
  
private:

  RefCountedPtr<mfis> m_mfis;

  Vector<EBLevelGrid> m_eblg;
};

#endif
