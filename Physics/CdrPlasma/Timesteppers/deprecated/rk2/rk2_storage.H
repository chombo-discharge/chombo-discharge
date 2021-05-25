/*!
  @file rk2_storage.H
  @brief Declaration of scratch storage for rk2
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _RK2_STORAGE_
#define _RK2_STORAGE_

#include "rk2.H"
#include "CD_AmrMesh.H"

#include <RefCountedPtr.H>

/*!
  @brief Class for handling extra data allocations that are necessary to advance a CdrSolver in an RK2 context. 
  @details This class needs storage for the intermediate state, k1, k2, and 3 irregular data holders for extrapolated fluxes,
  extrapolated densities, and extrapolates velocities. 
*/
class rk2::CdrStorage {
public:

  /*!
    @brief Weak constructor
  */
  CdrStorage();

  /*!
    @brief Constructor
  */
  CdrStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  ~CdrStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get phi
  */
  virtual EBAMRCellData& getPhi(){
    return m_phi;
  }

  /*!
    @brief Scratch state
  */
  virtual EBAMRCellData& getScratch(){
    return m_scratch;
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
  virtual EBAMRIVData& getEbState(){
    return m_scratchIV1;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& getEbVelo(){
    return m_scratchIV2;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& getEbFlux(){
    return m_scratchIV3;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& getEbGrad(){
    return m_scratchIV4;
  }

protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<AmrMesh> m_amr;

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
    @brief Extra storage required for transient RTE solvers
  */
  EBAMRCellData m_scratch;

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

  /*!
    @brief Scratch
  */
  EBAMRIVData m_scratchIV4;
};

/*!
  @brief Clas for handling extra storage requirements for advancing the Poisson solver with rk2. This only includes the 
  intermediate state
*/
class rk2::FieldStorage {
public:

  /*!
    @brief Constructor
  */
  FieldStorage();

  /*!
    @brief Full constructor
  */
  FieldStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  virtual ~FieldStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get intermediate state
  */
  virtual MFAMRCellData& getPhi(){
    return m_phi;
  }

  /*!
    @brief Get scratch storage
  */
  virtual MFAMRCellData& getScratch_phi(){
    return m_scratch_phi;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRCellData& getScratch_E(){
    return m_scratch_E;
  }

  /*!
    @brief Get the electric field storage
  */
  virtual EBAMRCellData& getElectricFieldCell(){
    return m_E_cell;
  }

  /*!
    @brief Get the electric field storage
  */
  virtual EBAMRFluxData& getElectricFieldFace(){
    return m_E_face;
  }

  /*!
    @brief Get the electric field storage
  */
  virtual EBAMRIVData& getElectricFieldEb(){
    return m_E_eb;
  }

protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<AmrMesh> m_amr;

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
    @brief Scratch phi. Only used for transient RTE solvers
  */
  MFAMRCellData m_scratch_phi;

  /*!
    @brief Scratch E. Only used for transient RTE solvers
  */
  EBAMRCellData m_scratch_E;

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
class rk2::RtStorage {
public:

  /*!
    @brief Weak constructor
  */
  RtStorage();

  /*!
    @brief Constructor
  */
  RtStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  ~RtStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();

  /*!
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get phi
  */
  virtual EBAMRCellData& getPhi(){
    return m_phi;
  }

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRIVData& getEbFlux(){
    return m_scratchIV;
  }

protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<AmrMesh> m_amr;

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

class rk2::SigmaStorage {
public:

  SigmaStorage();
  
  /*!
    @brief Constructor
  */
  SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);

  /*!
    @brief Destructor
  */
  ~SigmaStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_storage();
  
  /*!
    @brief Deallocate storage
  */
  virtual void deallocate_storage();

  /*!
    @brief Get phi
  */
  virtual EBAMRIVData& getPhi(){
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
    return m_k2;
  }
  
protected:

  /*!
    @brief AMR mesh
  */
  RefCountedPtr<AmrMesh> m_amr;

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
#include "CD_NamespaceFooter.H"
#endif
   
