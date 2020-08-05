/*!
  @file   godunov_storage.H
  @brief  Declaration of scratch storage for godunov
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _GODUNOV_STORAGE_
#define _GODUNOV_STORAGE_

#include "godunov.H"
#include "amr_mesh.H"

#include <RefCountedPtr.H>

namespace physics::cdr_plasma {
  
  class godunov::cdr_storage {
  public:

    cdr_storage();
    cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
		const std::string              a_realm,
		const phase::which_phase       a_phase,
		const int m_ncomp);
    
    ~cdr_storage();

    void allocate_storage();
    void deallocate_storage();

    EBAMRCellData& get_scratch()  { return m_scratch;    }
    EBAMRCellData& get_scratch2() { return m_scratch2;   }
    EBAMRCellData& get_scratch3() { return m_scratch3;   }
    EBAMRCellData& get_extrap()   { return m_cellExtr;  }
    EBAMRCellData& get_gradient() { return m_gradient;   }

    EBAMRIVData& get_eb_scratch1(){ return m_scratchIVs; }
    EBAMRIVData& get_eb_scratchD(){ return m_scratchIVD; }
    EBAMRIVData& get_eb_state(){    return m_scratchIV1; }
    EBAMRIVData& get_eb_velo(){     return m_scratchIV2; }
    EBAMRIVData& get_eb_flux(){     return m_scratchIV3; }
    EBAMRIVData& get_eb_grad(){     return m_scratchIV4; }

    EBAMRIFData& get_domain_state(){  return m_scratchIF1; }
    EBAMRIFData& get_domain_velo(){   return m_scratchIF2; }
    EBAMRIFData& get_domain_flux(){   return m_scratchIF3; }
    EBAMRIFData& get_domain_grad(){   return m_scratchIF4; }
  
  protected:
    std::string             m_realm;
    RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
    phase::which_phase m_phase;    // on this phase
    int m_ncomp;                   // with this many components

    EBAMRCellData m_scratch;  // Scratch data
    EBAMRCellData m_scratch2; // Scratch data
    EBAMRCellData m_scratch3; // Scratch data
    EBAMRCellData m_cellExtr; // Scratch data
    EBAMRCellData m_gradient; // Storage for holding gradient

    EBAMRIVData m_scratchIVs; // Scratch data
    EBAMRIVData m_scratchIVD; // Scratch data
    EBAMRIVData m_scratchIV1; // Scratch data
    EBAMRIVData m_scratchIV2; // Scratch data
    EBAMRIVData m_scratchIV3; // Scratch data
    EBAMRIVData m_scratchIV4; // Scratch data

    EBAMRIFData m_scratchIF1; // Scratch data
    EBAMRIFData m_scratchIF2; // Scratch data
    EBAMRIFData m_scratchIF3; // Scratch data
    EBAMRIFData m_scratchIF4; // Scratch data
  };

  class godunov::poisson_storage {
  public:

    poisson_storage();
    poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
		    const std::string              a_realm,
		    const phase::which_phase       a_phase,
		    const int m_ncomp);
    ~poisson_storage();

    void allocate_storage();
    void deallocate_storage();

    EBAMRCellData& get_E_cell()      {return m_E_cell;   }
    EBAMRIVData& get_E_eb()          {return m_E_eb;     }
    EBAMRIFData& get_E_domain()      {return m_E_dom;    }

  protected:
    std::string m_realm;
    RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
    phase::which_phase m_phase;    // on this phase
    int m_ncomp;                   // with this many components

    EBAMRCellData m_E_cell;        // Cell-centered E
    EBAMRIVData m_E_eb;            // EB-centered E
    EBAMRIFData m_E_dom;           // E-field on domain faces
  };

  class godunov::rte_storage {
  public:

    rte_storage();
    rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
		const std::string              a_realm,
		const phase::which_phase       a_phase,
		const int                      a_ncomp);
    ~rte_storage();

    void allocate_storage();
    void deallocate_storage();

    EBAMRIVData& get_eb_flux(){     return m_scratchIV; }
    EBAMRIFData& get_domain_flux(){ return m_scratchIF; }

  protected:

    std::string m_realm;
    RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
    phase::which_phase m_phase;    // on this phase
    int m_ncomp;                   // with this many components

    EBAMRIVData m_scratchIV;
    EBAMRIFData m_scratchIF;
  };

  class godunov::sigma_storage {
  public:

    sigma_storage();
    sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
		  const std::string              a_realm,
		  const phase::which_phase       a_phase,
		  const int                      a_ncomp);
    ~sigma_storage();
  
    void allocate_storage();
    void deallocate_storage();

    EBAMRIVData& get_scratch() {return m_scratch;}
  
  protected:

    std::string m_realm;
    RefCountedPtr<amr_mesh> m_amr; // Needed for allocation. Will allocate
    phase::which_phase m_phase;    // on this phase
    int m_ncomp;                   // with this many components

    EBAMRIVData m_scratch;
  };

}
#endif
