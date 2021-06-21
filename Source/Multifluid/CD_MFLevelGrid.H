/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFLevelGrid.H
  @brief  Declaration of a wrapper for wrapping multifluid EBLevelGrids
  @author Robert Marskar
  @todo   Remove old code when it's no longer needed (after multifluid_rework). 
*/

#ifndef CD_MFLevelGrid_H
#define CD_MFLevelGrid_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_Phases.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Wrapper class for holding multifluid EBLevelGrids
*/
class MFLevelGrid {
public:

  /*!
    @brief Weak constructor, must subsequently call define
  */
  MFLevelGrid();

  /*!
    @brief Full constructor
    @param[in] a_eblgs EBLevelGrids on both phases. 
  */
  MFLevelGrid(const std::map<Phase, EBLevelGrid>& a_eblgs);

  /*!
    @brief Define function
    @param[in] a_eblgs EBLevelGrids on both phases. 
  */
  void define(const std::map<Phase, EBLevelGrid>& a_eblgs);

  /*!
    @brief Get EBLevelGrid on specific phase
  */
  const EBLevelGrid& getEBLevelGrid(const Phase a_phase) const;

  /*!
    @brief Compute the interface region
    @param[in] a_box Grid box
    @param[in] a_dit Data index
  */
  IntVectSet getInterfaceRegion(const Box& a_box, const DataIndex& a_dit);

  /*!
    @brief Get number of phases
    @return Size of m_eblgs.
  */
  int getNumPhases() const;

  // Old stuff below here. 

  /*!
    @brief Full constructor. Does full allocation
  */
  MFLevelGrid(const DisjointBoxLayout&          a_dbl,
	      const ProblemDomain&              a_domain,
	      const int                         a_ebghost,
	      const RefCountedPtr<MultiFluidIndexSpace>&        a_multiFluidIndexSpace);

  /*!
    @brief Alternative full constructor - RefCounted storage so no allocation here
  */
  MFLevelGrid(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace, const Vector<EBLevelGrid>& a_eblg);

  /*!
    @brief Destructor
  */
  ~MFLevelGrid();

  /*!
    @brief Number of phases
  */
  virtual int numPhases() const;

  /*!
    @brief Define function
  */
  virtual void define(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace, const Vector<EBLevelGrid>& a_eblg);

  /*!
    @brief Set max refinement ratio
  */
  virtual void setMaxRefinementRatio(const int a_refRat);

  /*!
    @brief Get the index space
  */
  virtual const RefCountedPtr<MultiFluidIndexSpace>& getMfIndexSpace() const;

  /*!
    @brief Get the problem domain
  */
  virtual ProblemDomain getDomain() const;

  /*!
    @brief Get the grids
  */
  virtual DisjointBoxLayout getGrids() const;
    
  /*!
    @brief Get phase
  */
  virtual EBLevelGrid& getEBLevelGrid(int a_phase);

  /*!
    @brief Get phase
  */
  virtual const EBLevelGrid& getEBLevelGrid(int a_phase) const;

  /*!
    @brief Get interface region between two phases.
  */
  virtual IntVectSet interfaceRegion(const Box&       a_box,
				     const DataIndex& a_dit,
				     const int        a_phase1 = 0,
				     const int        a_phase2 = 1) const;

  /*!
    @brief Get the interface IntVect on the other side of the EB. Useful for when EB aligns with cell face
  */
  virtual bool interfacePair(IntVect& a_iv,
			     const IntVect&   a_iv_in,
			     const Box&       a_box,
			     const DataIndex& a_dit,
			     const int        a_phase1 = 0,
			     const int        a_phase2 = 1) const;
private:

  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  Vector<EBLevelGrid> m_eblg;

  /*!
    @brief EBLevelGrids on each phase
  */
  std::map<Phase, EBLevelGrid> m_eblgs;
};

#include <CD_NamespaceFooter.H>

#endif
