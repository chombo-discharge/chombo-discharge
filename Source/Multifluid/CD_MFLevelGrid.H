/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFLevelGrid.H
  @brief  Declaration of a wrapper for wrapping multifluid EBLevelGrids
  @author Robert Marskar
*/

#ifndef CD_MFLevelGrid_H
#define CD_MFLevelGrid_H

// Std includes
#include <map>

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Wrapper class for holding multifluid EBLevelGrids
*/
class MFLevelGrid
{
public:
  /*!
    @brief Weak constructor, must subsequently call define
  */
  MFLevelGrid();

  /*!
    @brief Full constructor. 
    @param[in] a_dbl     Grids
    @param[in] a_domain  Domains
    @param[in] a_ebghost Number of eb ghosts to fill
    @param[in] a_mfis    Multifluid index space
  */
  MFLevelGrid(const DisjointBoxLayout&                   a_dbl,
              const ProblemDomain&                       a_domain,
              const int                                  a_ebghost,
              const RefCountedPtr<MultiFluidIndexSpace>& a_mfis);

  /*!
    @brief Alternative full constructor - RefCounted storage so no allocation here
    @param[in] a_mfis  Multifluid index space
    @param[in] a_eblgs Levelgrids on each phase. 
  */
  MFLevelGrid(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace, const Vector<EBLevelGrid>& a_eblgs);

  /*!
    @brief Destructor
  */
  virtual ~MFLevelGrid();

  /*!
    @brief Get number of phases
  */
  virtual int
  numPhases() const;

  /*!
    @brief Define function
    @param[in] a_mfis  Multifluid index space
    @param[in] a_eblgs Levelgrids on each phase. 
  */
  virtual void
  define(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace, const Vector<EBLevelGrid>& a_eblgs);

  /*!
    @brief Set max refinement ratio
    @param[in] a_refRat Max refinement ratio. 
  */
  virtual void
  setMaxRefinementRatio(const int a_refRat);

  /*!
    @brief Get the index space
  */
  virtual const RefCountedPtr<MultiFluidIndexSpace>&
  getMfIndexSpace() const;

  /*!
    @brief Get the problem domain
  */
  virtual ProblemDomain
  getDomain() const;

  /*!
    @brief Get the grids
  */
  virtual DisjointBoxLayout
  getGrids() const;

  /*!
    @brief Get grids on specified phase
    @param[in] a_phase Phase
  */
  virtual EBLevelGrid&
  getEBLevelGrid(int a_phase);

  /*!
    @brief Get grids on specified phase
    @param[in] a_phase Phase
  */
  virtual const EBLevelGrid&
  getEBLevelGrid(int a_phase) const;

  /*!
    @brief Get interface region between two phases.
    @param[in] a_box    Input box
    @param[in] a_dit    Grid index
    @param[in] a_phase1 First phase
    @param[in] a_phase2 Second phase
    @return Returns the cells that are irregular on both phases. 
  */
  virtual IntVectSet
  interfaceRegion(const Box& a_box, const DataIndex& a_dit, const int a_phase1 = 0, const int a_phase2 = 1) const;

#if 0 // Taken out of chombo-discharge for now
  /*!
    @brief Get the interface IntVect on the other side of the EB. Useful for when EB aligns with cell face
  */
  virtual bool interfacePair(IntVect& a_iv,
			     const IntVect&   a_iv_in,
			     const Box&       a_box,
			     const DataIndex& a_dit,
			     const int        a_phase1 = 0,
			     const int        a_phase2 = 1) const;
#endif

private:
  /*!
    @brief Index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Grids on each level
  */
  Vector<EBLevelGrid> m_eblg;
};

#include <CD_NamespaceFooter.H>

#endif
