/*!
  @file adaptive_rkSSP_storage.H
  @brief Declaration of scratch storage for adaptive_rkSSP
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _ADAPTIVE_RKSSP_STORAGE_
#define _ADAPTIVE_RKSSP_STORAGE_

#include "adaptive_rkSSP.H"
#include "amr_mesh.H"

#include <RefCountedPtr.H>

class adaptive_rkSSP::cdr_storage {
public:

  cdr_storage();
  cdr_storage(const int m_stages, const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~cdr_storage();

  void allocate_storage();
  void deallocate_storage();

  EBAMRCellData& get_cache(){   return m_cache;      }
  EBAMRCellData& get_scratch(){ return m_scratch;    }
  EBAMRCellData& get_previous(){return m_previous;   }
  EBAMRCellData& get_error(){   return m_error;      }
  EBAMRCellData& get_gradient(){return m_gradient;   }

  EBAMRIVData& get_eb_state(){  return m_scratchIV1; }
  EBAMRIVData& get_eb_velo(){   return m_scratchIV2; }
  EBAMRIVData& get_eb_flux(){   return m_scratchIV3; }
  EBAMRIVData& get_eb_grad(){   return m_scratchIV4; }

  EBAMRIFData& get_domain_state(){  return m_scratchIF1; }
  EBAMRIFData& get_domain_velo(){   return m_scratchIF2; }
  EBAMRIFData& get_domain_flux(){   return m_scratchIF3; }
  EBAMRIFData& get_domain_grad(){   return m_scratchIF4; }
  
protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages

  EBAMRCellData m_cache;    // Cached solution
  EBAMRCellData m_scratch;  // Scratch data
  EBAMRCellData m_previous; // Previous solution (for each substep)
  EBAMRCellData m_error;    // Error estimator
  EBAMRCellData m_gradient; // Storage for holding gradient

  EBAMRIVData m_scratchIV1; // Scratch data
  EBAMRIVData m_scratchIV2; // Scratch data
  EBAMRIVData m_scratchIV3; // Scratch data
  EBAMRIVData m_scratchIV4; // Scratch data

  EBAMRIFData m_scratchIF1;
  EBAMRIFData m_scratchIF2;
  EBAMRIFData m_scratchIF3;
  EBAMRIFData m_scratchIF4;
};

class adaptive_rkSSP::poisson_storage {
public:

  poisson_storage();
  poisson_storage(const int m_stages, const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~poisson_storage();

  void allocate_storage();
  void deallocate_storage();
  
  MFAMRCellData& get_cache(){   return m_cache;   }
  EBAMRCellData& get_E_cell(){  return m_E_cell;  }
  EBAMRFluxData& get_E_face(){  return m_E_face;  }
  EBAMRIVData& get_E_eb(){      return m_E_eb;    }
  EBAMRIFData& get_E_domain(){  return m_E_dom;   }

protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages

  MFAMRCellData m_cache;         // Cached solution
  EBAMRCellData m_E_cell;        // Cell-centered E
  EBAMRFluxData m_E_face;        // Face-centered E
  EBAMRIVData m_E_eb;            // EB-centered E
  EBAMRIFData m_E_dom;           // E-field on domain faces
};

class adaptive_rkSSP::rte_storage {
public:

  rte_storage();
  rte_storage(const int a_stages, const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~rte_storage();

  void allocate_storage();
  void deallocate_storage();
  
  EBAMRCellData& get_cache(){ return m_cache; }
  EBAMRIVData& get_eb_flux(){ return m_scratchIV; }
  EBAMRIFData& get_domain_flux(){ return m_scratchIF; }

protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages


  EBAMRCellData m_cache;   // Cached solution
  EBAMRIVData m_scratchIV; // Used for EB-centered isotropic density
  EBAMRIFData m_scratchIF; // Used for domain-centered isotropic density
};

class adaptive_rkSSP::sigma_storage {
public:

  sigma_storage();
  sigma_storage(const int a_stages, const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~sigma_storage();
  
  void allocate_storage();
  void deallocate_storage();
  
  EBAMRIVData& get_cache(){   return m_cache; }
  EBAMRIVData& get_scratch(){ return m_scratch; }
  EBAMRIVData& get_previous(){return m_previous; }
  EBAMRIVData& get_error(){   return m_error; }

  
protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages

  EBAMRIVData m_cache;           // Cached solution
  EBAMRIVData m_scratch;         // Scratch data
  EBAMRIVData m_previous;        // Previous solution
  EBAMRIVData m_error;           // Error 
};

#endif
