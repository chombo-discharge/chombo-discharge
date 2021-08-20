/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxType.H
  @brief  Declaration of a class for identifying boxas as regular, cut, or covered. 
  @author Robert Marskar
*/

#ifndef CD_BoxType_H
#define CD_BoxType_H

// Chombo includes
#include <Box.H>
#include <LevelData.H>

#include <CD_NamespaceHeader.H>

/*!
  @brief Class for describing if a box is cut, regular, or covered. 
  @details This can be used in parallelas LevelData<BoxType>
*/
class BoxType {
public:

  /*!
    @brief Pre-allocatable or not
  */
  static int preAllocatable();

  /*!
    @brief Define box type
  */
  BoxType();

  /*!
    @brief Define box over a box
    @param[in] a_box Input box
    @param[in] comps Input components. 
    @details Redundant function object. Sets the box type to -1 => unknown box type
  */
  BoxType(const Box& box, int comps);

  /*!
    @brief Set box to covered (m_which = 0)
  */
  void setCovered();

  /*!
    @brief Set box to regular (m_which=1)
  */
  void setRegular();

  /*!
    @brief Set box to cutcell (m_which=2)
  */
  void setCutCell();

  /*!
    @brief Check if box is covered
    @return True if the box is covered
  */
  bool isCovered() const;

  /*!
    @brief Check if box is regular
    @return True if the box is regular
  */
  bool isRegular() const;

  /*!
    @brief Check if box is cut
    @return True if the box contains cut-cells
  */
  bool isCutCell() const;

  /*!
    @brief Define function
    @details Input parameters are ignored and the box type is set to m_which=-1
  */
  void define(const Box& box, int comps);

  /*!
    @brief Copy function. This sets this box type to be the input box's boxtype. 
    @param[in] R      Input region
    @param[in] Cd     Destination interval
    @param[in] source Source 
    @param[in] Cs     Source interval
    @details If the input/output boxes overlap we copy over the box type. Nothing else. 
  */
  void copy(const Box& R, const Interval& Cd, const BoxType& source, const Interval Cs);

  /*!
    @brief Copy function
    @param[in] a_regionFrom Region
    @param[in] a_Cdest      Interval
    @param[in] a_RegionTo   Box
    @param[in] a_src        Boxtype
    @param[in] a_Csrc       Interval
  */
  void copy(const Box& a_RegionFrom, const Interval& a_Cdest, const Box& a_RegionTo, const BoxType& a_src, const Interval& a_Csrc);

  /*!
    @brief Size function
    @param[in] R Box
    @param[in] comps Components
    @return size of an int (we are only communicated m_which)
  */  
  int size(const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
    @param[inout] buf   Buffer
    @param[in]    R     Input box
    @param[in]    comps Input components
    @note This write m_which to the buffer
  */
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
    @param[inout] buf   Buffer
    @param[in]    R     Input box
    @param[in]    comps Input components
    @note This sets m_which from buf
  */
  void linearIn(void* buf, const Box& R, const Interval& comps);

  /*!
    @brief Assignment function
    @param[in] a_src Source 
    @details This sets m_which to be a_src.m_which
  */
  void operator=(const BoxType& a_src);


  int m_which;
};

/*!
  @brief Factory class for BoxType
*/
class BoxTypeFactory : public DataFactory<BoxType> {
public:

  /*!
    @brief Default constructor
  */
  BoxTypeFactory();

  /*!
    @brief Destructor
  */  
  ~BoxTypeFactory();

  /*!
    @brief Factory method. Returns a new BoxType with m_which=-1
    @param[in] a_box    Input box
    @param[in] a_ncomps Number of compnents
    @param[in] a_dit    Grid idnex
  */
  virtual BoxType* create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const;
};

#include <CD_NamespaceFooter.H>

#endif
