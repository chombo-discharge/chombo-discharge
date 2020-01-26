/*!
  @file   BoxType.H
  @brief  Simple class for describing if a box is cut, regular, or covered
  @author Robert Marskar
*/

#ifndef _BOXTYPE_
#define _BOXTYPE_

#include "Box.H"
#include "LevelData.H"

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
  void copy(const Box& a_RegionFrom, const Interval& a_Cdest,  const Box& a_RegionTo, const BoxType& a_src, const Interval& a_Csrc);

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

protected:

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

#endif
