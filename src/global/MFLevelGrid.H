/*!
  @file MFLevelGrid
  @brief Declaration of a wrapper for wrapping multifluid EBLevelGrids
  @author Robert Marskar
  @date Dec. 2017
  @todo Segregate header and implementation files
*/

#ifndef _MFLEVELGRID_
#define _MFLEVELGRID_

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

  MFLevelGrid(const DisjointBoxLayout&          a_dbl,
	      const ProblemDomain&              a_domain,
	      const int                         a_ebghost,
	      const Vector<const EBIndexSpace*> a_ebis){
    m_eblg.resize(0);
    for (int i = 0; i < a_ebis.size(); i++){
      m_eblg.push_back(EBLevelGrid(a_dbl, a_domain, a_ebghost, a_ebis[i]));
    }
  }

  /*!
    @brief Constructor
  */
  MFLevelGrid(Vector<EBLevelGrid>& a_eblg){
    this->define(a_eblg);
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
  virtual void define(Vector<EBLevelGrid>& a_eblg){
    m_eblg = a_eblg;
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
    @brief Get EBISLayouts
  */
  virtual const Vector<EBISLayout> get_ebisl() const {
    MayDay::Abort("MFLevelGrid::get_ebisl() - code appers to be broken"); // Maybe we should return by reference
    Vector<EBISLayout> ebisl;
    for (int i = 0; i < m_eblg.size(); i++){
      ebisl.push_back(m_eblg[i].getEBISL());
    } // This appears to be broekN!
    return ebisl;
  }

  virtual const Vector<const EBIndexSpace*> get_ebis() const {
    Vector<const EBIndexSpace*> ebis;
    for (int i = 0; i < m_eblg.size(); i++){
      ebis.push_back(m_eblg[i].getEBIS());
    }

    return ebis;
  }
  
private:

  Vector<EBLevelGrid> m_eblg;
};

#endif
