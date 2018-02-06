/*!
  @file rk2_storage.H
  @brief Declaration of scratch storage for rk2
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _RK2_STORAGE_
#define _RK2_STORAGE_

#include "rk2.H"
#include "amr_mesh.H"

#include <RefCountedPtr.H>

/*!
  @brief Class for handling extra data allocations that are necessary to advance a cdr_solver in an RK2 context. 
  @details This class needs storage for the intermediate state, k1, k2, and 3 irregular data holders for extrapolated fluxes,
  extrapolated densities, and extrapolates velocities. 
*/
class rk2::cdr_storage {
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
    @brief Get phi
  */
  virtual EBAMRCellData& get_phi(){
    return m_phi;
  }

  /*!
    @brief Get k1
  */
  virtual EBAMRCellData& get_k1(){
    return m_k1;
  }

  /*!
    @brief Get k2
  */
  virtual EBAMRCellData& get_k2(){
    return m_k2;
  }

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
    @brief Intermediate state
  */
  EBAMRCellData m_phi;

  /*!
    @brief k1
  */
  EBAMRCellData m_k1;

  /*!
    @brief k2
  */
  EBAMRCellData m_k2;

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
  @brief Clas for handling extra storage requirements for advancing the Poisson solver with rk2. This only includes the 
  intermediate state
*/
class rk2::poisson_storage {
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
    @brief Get intermediate state
  */
  virtual MFAMRCellData& get_phi(){
    return m_phi;
  }

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
    @brief Intermediate state
  */
  MFAMRCellData m_phi;

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
  @brief Clas for handling extra storage requirements for advancing rte equations with rk2. This only includes the 
  intermediate state and the eb flux
*/
class rk2::rte_storage {
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
    @brief Get phi
  */
  virtual EBAMRCellData& get_phi(){
    return m_phi;
  }

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
    @brief k2
  */
  EBAMRCellData m_phi;

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV;
};

class rk2::sigma_storage {
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
    @brief Get phi
  */
  virtual EBAMRIVData& get_phi(){
    return m_phi;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_k1(){
    return m_k1;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_k2(){
    return m_k1;
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
    @brief State
  */
  EBAMRIVData m_phi;

  /*!
    @brief k1
  */
  EBAMRIVData m_k1;

  /*!
    @brief k2
  */
  EBAMRIVData m_k2;
};
#endif
   
