/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODESolver.H
  @brief  Encapsulation of an ODE solver on the mesh. 
  @author Robert Marskar
*/

#ifndef CD_MeshODESolver_H
#define CD_MeshODESolver_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for solving dy/dt = f on an AMR hierarchy. 
  @details The template parameter is the number of variables in y and f.
*/
template <size_t N = 1>
class MeshODESolver
{
public:
  /*!
    @brief Default constructor. Must subsequently set everything through public member functions.
  */
  MeshODESolver();

  /*!
    @brief Disallowed copy constructor
  */
  MeshODESolver(const MeshODESolver&) = delete;

  /*!
    @brief Disallowed move constructor
  */
  MeshODESolver(const MeshODESolver&&) = delete;

  /*!
    @brief Destructor
  */
  virtual ~MeshODESolver();

  /*!
    @brief Disallowed assignment operator
  */
  MeshODESolver&
  operator=(const MeshODESolver&) = delete;

  /*!
    @brief Disallowed move assignement operator
  */
  MeshODESolver&
  operator=(const MeshODESolver&&) = delete;

  /*!
    @brief Parse class options.
  */
  virtual void
  parseOptions() noexcept;

  /*!
    @brief Parse run-time configurable class options.
  */
  virtual void
  parseRuntimeOptions() noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  EBAMRCellData&
  getY() noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  const EBAMRCellData&
  getY() const noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  EBAMRCellData&
  getRHS() noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  const EBAMRCellData&
  getRHS() const noexcept;

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This copies m_y onto m_cache
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @details This linearly interpolates (potentially with limiters) m_y to the new grids. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

protected:
  /*!
    @brief Solver name
  */
  std::string m_name;

  /*!
    @brief Class name. 
    @details This will be different for different implementations of CdrSolver
  */
  std::string m_className;

  /*!
    @brief Solution vector (left hand side of equation). 
  */
  EBAMRCellData m_y;

  /*!
    @brief Source term (right-hand side of equation). 
  */
  EBAMRCellData m_rhs;

  /*!
    @brief Cache storage when doing regrids.
  */
  EBAMRCellData m_cache;

  /*!
    @brief Use slopes when regridding or not
  */
  bool m_regridSlopes;
};

#include <CD_NamespaceFooter.H>

#include <CD_MeshODESolver.H>

#endif
