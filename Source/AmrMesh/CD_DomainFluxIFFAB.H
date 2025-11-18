/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DomainFluxIFFAB.H
  @brief  Declaration of a BaseIFFAB wrapper that holds domain fluxes 
  @author Robert Marskar
*/

#ifndef CD_DomainFluxIFFAB_H
#define CD_DomainFluxIFFAB_H

// Chombo includes
#include <BaseIFFAB.H>
#include <ProblemDomain.H>
#include <EBISLayout.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding fluxes on domain boundaries
  @details The internal functionality of this class is to use two SpaceDim-sized tuples that hold the fluxes in each direction. Object in this class
  can be linearized and thus put in LevelData as LevelData<DomainFluxIFFAB>. 
  @note The BaseIFFABs are defined over the strip of cells abutting a domain boundary -- it it will also contain a layer of interior cells. 
*/
class DomainFluxIFFAB
{
public:
  static int
  preAllocatable()
  {
    return 1; // Symmetric messaging
  }

  /*!
    @brief Default constructor. Must subsequently call define
  */
  DomainFluxIFFAB();

  /*!
    @brief Full constructor. Calls the define function. 
    @param[in] a_domain  Problem domain
    @param[in] a_ebisbox EBISBox 
    @param[in] a_box     Cell-centered box in which we define the BaseIFFABs
    @param[in] a_nComp   Number of components
  */
  DomainFluxIFFAB(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_nComp);

  /*!
    @brief Destructor (does nothing)
  */
  ~DomainFluxIFFAB();

  /*!
    @brief Get the defining box
  */
  const Box&
  box() const;

  /*!
    @brief Get number of components
  */
  int
  nComp() const;

  /*!
    @brief Get the defining problem domain
  */
  const ProblemDomain&
  getDomain() const;

  /*!
    @brief Get the defining ebis box
  */
  const EBISBox&
  getEBISBox() const;

  /*!
    @brief Make a clone of an input data holder
    @param[in] a_copy Other domain flux
  */
  void
  define(const DomainFluxIFFAB& a_copy);

  /*!
    @brief Full define function
    @param[in] a_domain  Problem domain
    @param[in] a_ebisbox EBISBox 
    @param[in] a_box     Cell-centered box in which we define the BaseIFFABs
    @param[in] a_nComp   Number of components
  */
  void
  define(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_nComp);

  /*!
    @brief Get the BaseIFFAB on a particular domain side
    @param[in] a_dir  Coordinate direction
    @param[in] a_side Side (lo/hi)
  */
  BaseIFFAB<Real>&
  operator()(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Get the BaseIFFAB on a particular domain side
    @param[in] a_dir  Coordinate direction
    @param[in] a_side Side (lo/hi)
  */
  const BaseIFFAB<Real>&
  operator()(const int a_dir, const Side::LoHiSide a_side) const;

  /*!
    @brief Linearization function.
    @return Returns the sum of the BaseIFFAB sizes
  */
  int
  size(const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function.
    @details This reuses BaseIFFAB linearization function
  */
  void
  linearOut(void* buf, const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function.
    @details This reuses BaseIFFAB linearization function
  */
  void
  linearIn(void* buf, const Box& R, const Interval& comps);

  /*! 
    @brief Copy functiopn
    @details Copies from a subsection of one box into another.
    Assumes the boxes are both in the same index space, and that
    box R is completely contained in both the src and destination
    boxes.  This function equired by BoxLayoutData 
  */
  void
  copy(const Box& Rfrom, const Interval& Cdest, const Box& Rto, const DomainFluxIFFAB& src, const Interval& Csrc);

protected:
  /*!
    @brief Defining domain
  */
  ProblemDomain m_domain;

  /*!
    @brief Defining box
  */
  Box m_box;

  /*!
    @brief EBISBox
  */
  EBISBox m_ebisbox;

  /*!
    @brief Number of components
  */
  int m_nComp;

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief BaseIFFABs on Side::Lo sides
  */
  Tuple<BaseIFFAB<Real>*, CH_SPACEDIM> m_fluxLo;

  /*!
    @brief BaseIFFABs on Side::Hi sides
  */
  Tuple<BaseIFFAB<Real>*, CH_SPACEDIM> m_fluxHi;

  /*!
    @brief Set default values
  */
  void
  setDefaultValues();

  /*!
    @brief Clear object
  */
  void
  clear();
};

#include <CD_NamespaceFooter.H>

#endif
