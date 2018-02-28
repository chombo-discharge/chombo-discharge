/*!
  @file euler_f_storage.H
  @brief Declaration of scratch storage for euler_f
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _EULER_F_STORAGE_
#define _EULER_F_STORAGE_

#include "euler_f.H"
#include "amr_mesh.H"

#include <RefCountedPtr.H>

/*!
  @brief Class for handling extra data allocations that are necessary to advance a cdr_solver for euler_f
*/
class euler_f::cdr_storage {
public:

  /*!
    @brief Weak constructor
  */
  cdr_storage();

  /*!
    @brief Constructor
  */
  cdr_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  ~cdr_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_eb_state(){
    return m_scratchIV1;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_eb_velo(){
    return m_scratchIV2;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_eb_flux(){
    return m_scratchIV3;
  }

protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Number of componenents
  */
  int m_ncomp;

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV1;

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV2;

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV3;
};

/*!
  @brief Class for handling extra storage requirements for advancing the Poisson solver with euler_f. 
*/
class euler_f::poisson_storage {
public:

  /*!
    @brief Constructor
  */
  poisson_storage();

  /*!
    @brief Full constructor
  */
  poisson_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  virtual ~poisson_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Get the electric field storage
  */
  virtual EBAMRCellData& get_E_cell(){
    return m_E_cell;
  }

  /*!
    @brief Get the electric field storage
  */
  virtual EBAMRFluxData& get_E_face(){
    return m_E_face;
  }

    /*!
    @brief Get the electric field storage
  */
  virtual EBAMRIVData& get_E_eb(){
    return m_E_eb;
  }

protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Number of componenents
  */
  int m_ncomp;

  /*!
    @brief Intermediate field
  */
  EBAMRCellData m_E_cell;

  /*!
    @brief Field on face centers
  */
  EBAMRFluxData m_E_face;

  /*!
    @brief 
  */
  EBAMRIVData m_E_eb;
};

/*!
  @brief Class for handling extra storage requirements for advancing rte equations with euler_f.
*/
class euler_f::rte_storage {
public:

  /*!
    @brief Weak constructor
  */
  rte_storage();

  /*!
    @brief Constructor
  */
  rte_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  ~rte_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_eb_flux(){
    return m_scratchIV;
  }

protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Number of componenents
  */
  int m_ncomp;

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV;
};

/*!
  @brief Class for ahndling extra storage requirements for advancing the sigma equations with euler_f
*/
class euler_f::sigma_storage {
public:

  sigma_storage();
  
  /*!
    @brief Constructor
  */
  sigma_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  ~sigma_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_rhs(){
    return m_rhs;
  }
  
protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Number of componenents
  */
  int m_ncomp;

  /*!
    @brief Right hand side
  */
  EBAMRIVData m_rhs;
};
#endif
