/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoxType.H
  @brief  Declaration of a class for identifying boxas as regular, cut, or covered. Used with ScanSHop.
  @author Robert Marskar
*/

#ifndef CD_BoxType_H
#define CD_BoxType_H

// Chombo includes
#include <Box.H>
#include <LevelData.H>

#include <CD_NamespaceHeader.H>

/*!
  @brief Class for describing if a box is cut, regular, or covered
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
  */
  BoxType(const Box& box, int comps);

  /*!
    @brief Set box to covered
  */
  void setCovered();

  /*!
    @brief Set box to regular
  */
  void setRegular();

  /*!
    @brief Set box to cutcell
  */
  void setCutCell();

  /*!
    @brief Check if box is covered
  */
  bool isCovered() const;

  /*!
    @brief Check if box is regular
  */
  bool isRegular() const;

  /*!
    @brief Check if box is cut
  */
  bool isCutCell() const;

  /*!
    @brief Define function
  */
  void define(const Box& box, int comps);

  /*!
    @brief Copy function
  */
  void copy(const Box& R, const Interval& Cd, const BoxType& source, const Interval Cs);

  /*!
    @brief Copy function
  */
  void copy(const Box& a_RegionFrom, const Interval& a_Cdest, const Box& a_RegionTo, const BoxType& a_src, const Interval& a_Csrc);

  /*!
    @brief Size function
  */  
  int size(const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
  */
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  /*!
    @brief Linearization function
  */
  void linearIn(void* buf, const Box& R, const Interval& comps);

  /*!
    @brief Copy function
  */
  void operator=(const BoxType& a_src);

  //protected:

  int m_which;
};

/*!
  @brief Factory class for BoxType
*/
class BoxTypeFactory : public DataFactory<BoxType> {
public:

  BoxTypeFactory();

  ~BoxTypeFactory();

  virtual BoxType* create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const;
};

#include <CD_NamespaceFooter.H>

#endif
