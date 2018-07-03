/*!
  @file   implicit_trapezoidal_storage.H
  @brief  Declaration of scratch storage for implicit_trapezoidal
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _IMPLICIT_trapezoidal_STORAGE_
#define _IMPLICIT_trapezoidal_STORAGE_

#include "implicit_trapezoidal.H"
#include "amr_mesh.H"

#include <RefCountedPtr.H>

/*!
  @brief Storage class for evolving cdr equations in an implicit trapezoidal context
*/
class implicit_trapezoidal::cdr_storage {
public:

  /*!
    @brief Constructor
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
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get cached state
  */
  virtual EBAMRCellData& get_cache(){
    return m_cache;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRCellData& get_scratch1(){
    return m_scratch1;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRCellData& get_scratch2(){
    return m_scratch2;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRCellData& get_scratch3(){
    return m_scratch3;
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

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& get_eb_grad(){
    return m_scratchIV4;
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
    @brief Cached state
  */
  EBAMRCellData m_cache;

  /*!
    @brief Scratch
  */
  EBAMRCellData m_scratch1;

  /*!
    @brief Scratch
  */
  EBAMRCellData m_scratch2;

  /*!
    @brief Scratch
  */
  EBAMRCellData m_scratch3;

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

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV4;
};

/*!
  @brief Storage class for evolving rte equations in an implicit trapezoidal context
*/
class implicit_trapezoidal::rte_storage {
public:

  /*!
    @brief Constructor
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
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get cached state
  */
  virtual EBAMRCellData& get_cache(){
    return m_cache;
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
    @brief Cached state
  */
  EBAMRCellData m_cache;

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV;
};

/*!
  @brief Storage class for evolving the Poisson equation in an implicit trapezoidal context
*/
class implicit_trapezoidal::poisson_storage {
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
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get cached state
  */
  virtual MFAMRCellData& get_cache(){
    return m_cache;
  }

  /*!
    @brief Get cell-centered field storage
  */
  virtual EBAMRCellData& get_E_cell(){
    return m_E_cell;
  }

  /*!
    @brief Get face-centered field storage
  */
  virtual EBAMRFluxData& get_E_face(){
    return m_E_face;
  }

  /*!
    @brief Get EB-centered field storage
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
    @brief Cached state
  */
  MFAMRCellData m_cache;

  /*!
    @brief Intermediate field
  */
  EBAMRCellData m_E_cell;

  /*!
    @brief Field on face centers
  */
  EBAMRFluxData m_E_face;

  /*!
    @brief Field on the EB
  */
  EBAMRIVData m_E_eb;
};

/*!
  @brief Storage class for evolving the sigma equation in an implicit trapezoidal context
*/
class implicit_trapezoidal::sigma_storage {
public:

  /*!
    @brief Constructor
  */
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
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get cached state
  */
  virtual EBAMRIVData& get_cache(){
    return m_cache;
  }

  /*!
    @brief Get scratch
  */
  virtual EBAMRIVData& get_scratch1(){
    return m_scratch1;
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
    @brief Cached solution
  */
  EBAMRIVData m_cache;

  /*!
    @brief Scratch storage
  */
  EBAMRIVData m_scratch1;
  
};
#endif
