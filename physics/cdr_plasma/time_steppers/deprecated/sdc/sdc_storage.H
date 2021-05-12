/*!
  @file   sdc_storage.H
  @brief  Declaration of scratch storage for sdc
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _SDC_STORAGE_
#define _SDC_STORAGE_

#include "sdc.H"
#include "CD_AmrMesh.H"

#include <RefCountedPtr.H>

class sdc::cdr_storage {
public:

  cdr_storage();
  cdr_storage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~cdr_storage();

  void allocate_storage(const int a_p);
  void deallocate_storage();

  EBAMRCellData& get_scratch()  { return m_scratch;    }
  EBAMRCellData& get_scratch2() { return m_scratch2;   }
  EBAMRCellData& get_error()    { return m_error;      }
  EBAMRCellData& get_gradient() { return m_gradient;   }
  EBAMRCellData& get_old()      { return m_old;        }
  EBAMRCellData& get_scratchD() { return m_scratchD;   }
  EBAMRCellData& get_divD()     { return m_divD;       }
  EBAMRCellData& get_divF()     { return m_divF;       }
  EBAMRCellData& get_cell_mob() { return m_cell_mob;   }
  EBAMRFluxData& get_face_mob() { return m_face_mob;   }
  
  Vector<EBAMRCellData>& get_phi()     { return m_phi; }
  Vector<EBAMRCellData>& get_FA()      { return m_FA;  }
  Vector<EBAMRCellData>& get_FD()      { return m_FD;  }
  Vector<EBAMRCellData>& get_FR()      { return m_FR;  }
  Vector<EBAMRCellData>& get_F()       { return m_F;   }

  EBAMRIVData& get_eb_mob()     { return m_eb_mob;     }
  EBAMRIVData& get_eb_scratch1(){ return m_scratchIVs; }
  EBAMRIVData& get_eb_scratchD(){ return m_scratchIVD; }
  EBAMRIVData& get_eb_state()   { return m_scratchIV1; }
  EBAMRIVData& get_eb_velo()    { return m_scratchIV2; }
  EBAMRIVData& get_eb_flux()    { return m_scratchIV3; }
  EBAMRIVData& get_eb_grad()    { return m_scratchIV4; }

  EBAMRIFData& get_domain_state(){  return m_scratchIF1; }
  EBAMRIFData& get_domain_velo(){   return m_scratchIF2; }
  EBAMRIFData& get_domain_flux(){   return m_scratchIF3; }
  EBAMRIFData& get_domain_grad(){   return m_scratchIF4; }
  
protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_p;                       // SDC order
  bool m_subcycle;               // Extra storage for subcycling

  EBAMRCellData m_scratch;  // Scratch data
  EBAMRCellData m_scratch2; // Scratch data
  EBAMRCellData m_error;    // Error estimator
  EBAMRCellData m_gradient; // Storage for holding gradient
  EBAMRCellData m_old;      // "Old solution". Only used for subcycling
  EBAMRCellData m_scratchD; // SpaceDim scratch storage
  EBAMRCellData m_divF;     // Storage for DivF
  EBAMRCellData m_divD;     // Storage for DivD
  EBAMRCellData m_cell_mob; // Cell-centered mobility
  EBAMRFluxData m_face_mob; // Face-centered mobility
  
  Vector<EBAMRCellData> m_phi;      // phi^k in the  SDC notation
  Vector<EBAMRCellData> m_FA;       // FA in the SDC notation
  Vector<EBAMRCellData> m_FD;       // FD in the SDC notation
  Vector<EBAMRCellData> m_FR;       // FR in the SDC notation
  Vector<EBAMRCellData> m_F;        // F = FA + FD + FR

  EBAMRIVData m_eb_mob;     // EB mobility
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

class sdc::poisson_storage {
public:

  poisson_storage();
  poisson_storage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~poisson_storage();

  void allocate_storage(const int a_p);
  void deallocate_storage();

  MFAMRCellData& get_previous()    {return m_previous; }
  EBAMRCellData& get_E_cell()      {return m_E_cell;   }
  EBAMRFluxData& get_E_face()      {return m_E_face;   }
  EBAMRIVData& get_E_eb()          {return m_E_eb;     }
  EBAMRIFData& get_E_domain()      {return m_E_dom;    }

protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_p;                   // SDC order

  MFAMRCellData m_previous;      // Solution at previous time step
  EBAMRCellData m_E_cell;        // Cell-centered E
  EBAMRFluxData m_E_face;        // Face-centered E
  EBAMRIVData m_E_eb;            // EB-centered E
  EBAMRIFData m_E_dom;           // E-field on domain faces
};

class sdc::rte_storage {
public:

  rte_storage();
  rte_storage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~rte_storage();

  void allocate_storage(const int a_p);
  void deallocate_storage();

  Vector<EBAMRCellData>& get_phi()   {return m_phi;       }
  EBAMRCellData& get_previous()      {return m_previous;  }
  EBAMRIVData& get_eb_flux()         {return m_scratchIV; }
  EBAMRIFData& get_domain_flux()     {return m_scratchIF; }

protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_p;                   // SDC order

  Vector<EBAMRCellData> m_phi;
  EBAMRCellData m_previous;  // Solution at previous time step
  EBAMRIVData m_scratchIV;   // Used for EB-centered isotropic density
  EBAMRIFData m_scratchIF;   // Used for domain-centered isotropic density
};

class sdc::sigma_storage {
public:

  sigma_storage();
  sigma_storage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~sigma_storage();
  
  void allocate_storage(const int a_p);
  void deallocate_storage();

  EBAMRIVData& get_scratch(){ return m_scratch; }
  EBAMRIVData& get_error(){   return m_error; }

  Vector<EBAMRIVData>& get_sigma(){ return m_sigma;}
  Vector<EBAMRIVData>& get_Fnew(){  return m_Fnew; }
  Vector<EBAMRIVData>& get_Fold(){  return m_Fold; }

  
protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_p;                   // SDC order

  EBAMRIVData m_scratch;         // Scratch data
  EBAMRIVData m_error;           // Error

  Vector<EBAMRIVData> m_sigma;   // sigma^k in the SDC notation
  Vector<EBAMRIVData> m_Fnew;    // New operator lsopes
  Vector<EBAMRIVData> m_Fold;    // Lagged operator slopes
#include "CD_NamespaceFooter.H"

#endif
