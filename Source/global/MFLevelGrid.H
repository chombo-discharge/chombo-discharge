/*!
  @file   MFLevelGrid.H
  @brief  Declaration of a wrapper for wrapping multifluid EBLevelGrids
  @author Robert Marskar
  @date   Dec. 2017
  @todo   Segregate header and implementation files
*/

#ifndef _MFLEVELGRID_
#define _MFLEVELGRID_

#include "mfis.H"

#include "EBLevelGrid.H"

#include "CD_NamespaceHeader.H"

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
    @brief Full constructor. Does full allocation
  */
  MFLevelGrid(const DisjointBoxLayout&          a_dbl,
	      const ProblemDomain&              a_domain,
	      const int                         a_ebghost,
	      const RefCountedPtr<mfis>&        a_mfis);

  /*!
    @brief Alternative full constructor - RefCounted storage so no allocation here
  */
  MFLevelGrid(const RefCountedPtr<mfis>& a_mfis, const Vector<EBLevelGrid>& a_eblg);

  /*!
    @brief Destructor
  */
  ~MFLevelGrid();

  /*!
    @brief Number of phases
  */
  virtual int num_phases() const;

  /*!
    @brief Define function
  */
  virtual void define(const RefCountedPtr<mfis>& a_mfis, const Vector<EBLevelGrid>& a_eblg);

  /*!
    @brief Get the index space
  */
  virtual const RefCountedPtr<mfis>& get_mfis() const;

  /*!
    @brief Get the problem domain
  */
  virtual ProblemDomain get_domain() const;

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
  virtual IntVectSet interface_region(const Box&       a_box,
				      const DataIndex& a_dit,
				      const int        a_phase1 = 0,
				      const int        a_phase2 = 1) const;

  /*!
    @brief Get the interface IntVect on the other side of the EB. Useful for when EB aligns with cell face
  */
  virtual bool interface_pair(IntVect& a_iv,
			      const IntVect&   a_iv_in,
			      const Box&       a_box,
			      const DataIndex& a_dit,
			      const int        a_phase1 = 0,
			      const int        a_phase2 = 1) const;
private:

  RefCountedPtr<mfis> m_multifluidIndexSpace;

  Vector<EBLevelGrid> m_eblg;
};
#include "CD_NamespaceFooter.H"

#endif
