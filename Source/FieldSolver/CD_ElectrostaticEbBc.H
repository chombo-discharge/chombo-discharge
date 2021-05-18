/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ElectrostaticEbBc.H
  @brief  Data holder for holding the voltage at electrodes.
  @author Robert Marskar
*/

#ifndef CD_ElectrostaticEbBc_H
#define CD_ElectrostaticEbBc_H

// Std includes
#include <functional>
#include <vector>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <electrode.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which maps boundary condition types to a side and direction. 
  @details This class is meant to hold boundary conditions for use in various FieldSolver code. It provides no functionality otherwise,
  and has no dependencies, other than standard C++14. 
*/
class ElectrostaticEbBc {
public:

  /*!
    @brief Function which maps f(R^3,t) : R. Used for setting the potential at an electrode. 
  */
  using BcFunction = std::function<Real(const RealVect a_position, const Real a_time)>;
  
  /*!
    @brief Initializing constructor. This does nothing. 
  */
  ElectrostaticEbBc();

  /*!
    @brief Destructor. This does nothing. 
  */
  virtual ~ElectrostaticEbBc();

  /*!
    @brief Clears the object. 
  */
  void clear();

  /*!
    @brief Add a boundary condition
    @param[in] a_electrode  The electrode.
    @param[in] a_bcFunction Boundary condition function.
    @details This appends an EbBc to the end of m_bcFunctions. 
    @note If the boundary condition already exists for a_electrode, it will be overwritten. 
  */
  void addEbBc(const electrode& a_electrode, const BcFunction& a_bcFunction);

  /*!
    @brief Set the potential for a_electrode to be given by a_bcFunction
    @param[in] a_electrode The electrode index in the vector. 
    @param[in] a_bcFunction Boundary condition function.
  */
  void setEbBc(const int a_electrode, const BcFunction& a_bcFunction);

  /*!
    @brief Return boundary condition associated with a particular electrode.
    @param[in] a_electrode Electrode
    @return Boundary condition associated with a_electrode. 
  */
  BcFunction& getBc(const int a_electrode);

  /*!
    @brief Return boundary condition associated with a particular electrode.
    @param[in] a_electrode Electrode
    @return Boundary condition associated with a_electrode.
  */
  const BcFunction& getBc(const int a_electrode) const;

  /*!
    @brief Get all boundary conditions
    @return Returns modifiable version of m_bcFunctions
  */
  std::vector<std::pair<electrode, BcFunction> >& getBcs();

    /*!
    @brief Get all boundary conditions
    @return Returns non-modifiable version of m_bcFunctions
  */
  const std::vector<std::pair<electrode, BcFunction> >& getBcs() const;

protected:

  /*!
    @brief Data structure which holds electrode and associated potential functions. 
  */
  std::vector<std::pair<electrode, BcFunction> > m_bcFunctions;
};

#include <CD_NamespaceFooter.H>

#endif
