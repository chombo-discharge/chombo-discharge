/*!
  @file multirate_eulerf_storage.H
  @brief Declaration of scratch storage for multirate_eulerf
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _MULTIRATE_EULERF_STORAGE_
#define _MULTIRATE_EULERF_STORAGE_

#include "multirate_eulerf.H"
#include "amr_mesh.H"

#include <RefCountedPtr.H>

class multirate_eulerf::cdr_storage {
public:

  cdr_storage();
  cdr_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~cdr_storage();

  void allocate_storage();
  void deallocate_storage();

  EBAMRCellData& get_cache(){   return m_cache;      }
  EBAMRCellData& get_scratch(){ return m_scratch;    }

  EBAMRIVData& get_eb_state(){  return m_scratchIV1; }
  EBAMRIVData& get_eb_velo(){   return m_scratchIV2; }
  EBAMRIVData& get_eb_flux(){   return m_scratchIV3; }
  EBAMRIVData& get_eb_grad(){   return m_scratchIV4; }

protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components

  EBAMRCellData m_cache;    // Cached solution
  EBAMRCellData m_scratch;  // Scratch data

  EBAMRIVData m_scratchIV1; // Scratch data
  EBAMRIVData m_scratchIV2; // Scratch data
  EBAMRIVData m_scratchIV3; // Scratch data
  EBAMRIVData m_scratchIV4; // Scratch data
};

class multirate_eulerf::poisson_storage {
public:

  poisson_storage();
  poisson_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~poisson_storage();

  void allocate_storage();
  void deallocate_storage();
  
  MFAMRCellData& get_cache(){  return m_cache;  }
  EBAMRCellData& get_E_cell(){ return m_E_cell; }
  EBAMRFluxData& get_E_face(){ return m_E_face; }
  EBAMRIVData& get_E_eb(){     return m_E_eb;   }

protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components

  MFAMRCellData m_cache;         // Cached solution
  EBAMRCellData m_E_cell;        // Cell-centered E
  EBAMRFluxData m_E_face;        // Face-centered E
  EBAMRIVData m_E_eb;            // EB-centered E
};

class multirate_eulerf::rte_storage {
public:

  rte_storage();
  rte_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~rte_storage();

  void allocate_storage();
  void deallocate_storage();
  
  EBAMRCellData& get_cache(){ return m_cache; }
  EBAMRIVData& get_eb_flux(){ return m_scratchIV; }

protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components


  EBAMRCellData m_cache;   // Cached solution
  EBAMRIVData m_scratchIV; // Used for EB-centered isotropic density
};

class multirate_eulerf::sigma_storage {
public:

  sigma_storage();
  sigma_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~sigma_storage();
  
   void allocate_storage();
   void deallocate_storage();

   EBAMRIVData& get_cache(){   return m_cache; }
   EBAMRIVData& get_scratch(){ return m_scratch; }

  
protected:

  RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  

  EBAMRIVData m_cache;           // Cached solution
  EBAMRIVData m_scratch;         // Scratch data
};

#endif
