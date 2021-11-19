/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbGhostCloud.H
  @brief  A class for handling particle deposition clouds near refinement boundaries. 
  @author Robert Marskar
*/

#ifndef CD_EbGhostCloud_H
#define CD_EbGhostCloud_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBLevelGrid.H>
#include <EBCellFAB.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief A class for handling particle deposition clouds that hang over refinement boundaries. This class does not include any deposition functionality but provides
  buffers for moving mass between two levels. 
  @details To use this class, the user can directly call addFineGhostsToCoarse to add the contribution of the fine level's particles to the coarse level. This
  will add the part of the particle clouds that hang over the refinement boundary to the coarse level. Likewise, the user can get a buffer (getFiCoBuffer) for
  getting a temporary data holder on which he can use special deposition procedures.
  @todo This class currently does not work with multi-valued cells (like most of our particle code). Moreover, we actually refine the entire coarse level so that
  users can deposit particles without moving them around. This allocates a lot of mesh data and we might want to reconsider the design of this class. Indeed, there
  are other ways of doing this kind of deposition. 
*/
class EbGhostCloud {
public:

  /*!
    @brief Default constructor. Leaves object in undefined state and subsequently requires a call to define. 
  */
  EbGhostCloud();

  /*!
    @brief Full constructor. Calls define. 
    @param[in] a_eblgCoar   Coarse grids
    @param[in] a_eblgFine   Fine grids
    @param[in] a_refRat     Refinement factory between levels
    @param[in] a_nComp      Number of components allocated in buffer data holders. 
    @param[in] a_ghost      Number of ghost cells in mesh data holders. 
  */
  EbGhostCloud(const EBLevelGrid& a_eblgCoar,
	       const EBLevelGrid& a_eblgFine,
	       const int          a_refRat,
	       const int          a_nComp,
	       const int          a_ghost);

  /*!
    @brief Destructor (does nothing)
  */
  ~EbGhostCloud();

  /*!
    @brief Define function
    @param[in] a_eblgCoar   Coarse grids
    @param[in] a_eblgFine   Fine grids
    @param[in] a_refRat     Refinement factory between levels
    @param[in] a_nComp      Number of components allocated in buffer data holders. 
    @param[in] a_ghost      Number of ghost cells in mesh data holders. 
  */
  void define(const EBLevelGrid& a_eblgCoar,
	      const EBLevelGrid& a_eblgFine,
	      const int          a_refRat,
	      const int          a_nComp,
	      const int          a_ghost);

  /*!
    @brief Coarsen the fine-level ghost cells and add them to the coarse level. 
    @details On output, the contests in the ghost cells on the refinement boundary are added from the fine level to the coarse level. This routine
    is used when the fine-level particles have particle clouds that deposit into the ghost cells outside the fine level. This routine coarsens that
    data and adds the mass to the coarse level. 
    @param[inout] a_coarData Coarse mesh data. 
    @param[in]    a_fineData Fine mesh data. 
  */
  void addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData) const;

  /*!
    @brief Add the ghost cells on the refined coarse level to the fine level. 
    @details During this procedure the user will have deposited the coarse-level particles onto a_fiCoData, which is the refined coarse data holder. This routine
    will add the contents in that data holder onto the fine level. 
    @param[inout] a_fineData Fine mesh data
    @param[in]    a_fiCoData Refined coarse data buffer. Users will have deposited the mesh data onto this 
  */
  void addFiCoDataToFine(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_fiCoData) const;

  /*!
    @brief Get the refined coarse buffer
  */
  LevelData<EBCellFAB>& getFiCoBuffer() const;

  /*!
    @brief Get the refined version of the coarse eblevelgrid
  */
  const EBLevelGrid& getEblgFiCo() const;

protected:

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Coarse EBLevelGrid
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Fine EBLevelGrid
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Refined coarse EBLevelGrid
  */
  EBLevelGrid m_eblgFiCo;

  /*!
    @brief Coarsening of m_eblgFine
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Copier for copying data from the refined coarse level to the fine level.
  */
  Copier m_copierFiCoToFine;

  /*!
    @brief Copier for copying data from the coarsened fine level to the coarse level.
  */
  Copier m_copierCoFiToCoar;

  /*!
    @brief Data defined on the coarsened fine grid
  */
  mutable LevelData<EBCellFAB> m_bufferCoFi;

  /*!
    @brief Data defined on the refined coarse grid
  */
  mutable LevelData<EBCellFAB> m_bufferFiCo;

  /*!
    @brief Scratch storage for fine grid data
  */
  mutable LevelData<EBCellFAB> m_bufferFine;

  /*!
    @brief VoFIterator for fine-grid irregular ghost cells. 
  */
  mutable LayoutData<VoFIterator> m_vofIterFineGhosts;

  /*!
    @brief VoFIterator for coarse grid cells that are coarsenings of fine-grid ghost irregular cells.
    @details Essentially a coarsening of m_vofIterFineGhosts.
  */
  mutable LayoutData<VoFIterator> m_vofIterCoFiGhosts;  

  /*!
    @brief Refinement ratio between the two levels
  */
  int m_refRat;

  /*!
    @brief Maximum number of components
  */
  int m_nComp;

  /*!
    @brief Required number of ghost cells
  */
  int m_ghost;

  /*!
    @brief Define the vof iterators
  */
  void defineVoFIterators();
};

#include <CD_NamespaceFooter.H>

#endif
