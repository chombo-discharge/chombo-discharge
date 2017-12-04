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
  const EBLevelGrid& get_eblg(int a_phase) const {
    return m_eblg[a_phase];
  }

  /*!
    @brief Get EBISLayouts
  */
  Vector<EBISLayout> get_ebisl() {
    Vector<EBISLayout> ebisl;
    for (int i = 0; i < m_eblg.size(); i++){
      ebisl.push_back(m_eblg[i].getEBISL());
    }
    return ebisl;
  }
  
private:

  Vector<EBLevelGrid> m_eblg;
};

#endif
