/*!
  @file   strang2_storage.H
  @brief  Declaration of scratch storage for strang2
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _STRANG2_STORAGE_
#define _STRANG2_STORAGE_

#include "strang2.H"
#include "CD_AmrMesh.H"

#include <RefCountedPtr.H>

class strang2::CdrStorage {
public:

  CdrStorage();
  CdrStorage(const int m_stages, const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~CdrStorage();

  void allocate_storage();
  void deallocate_storage();
  void allocate_extra_storage(const int a_num_extra);
  void deallocate_extra_storage();

  EBAMRCellData& get_cache(){   return m_cache;      }
  EBAMRCellData& get_backup(){  return m_backup;     }
  EBAMRCellData& getScratch(){ return m_scratch;    }
  EBAMRCellData& getPrevious(){return m_previous;   }
  EBAMRCellData& getError(){   return m_error;      }
  EBAMRCellData& getGradient(){return m_gradient;   }
  
  Vector<EBAMRCellData*> get_extra_storage(){return m_extra_storage;}

  EBAMRIVData& getEbState(){  return m_scratchIV1; }
  EBAMRIVData& getEbVelo(){   return m_scratchIV2; }
  EBAMRIVData& getEbFlux(){   return m_scratchIV3; }
  EBAMRIVData& getEbGrad(){   return m_scratchIV4; }

  EBAMRIFData& getDomainState(){  return m_scratchIF1; }
  EBAMRIFData& getDomainVelo(){   return m_scratchIF2; }
  EBAMRIFData& getDomainFlux(){   return m_scratchIF3; }
  EBAMRIFData& getDomainGrad(){   return m_scratchIF4; }
  
protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages
  bool m_has_extra;              // 

  EBAMRCellData m_cache;    // Cached solution
  EBAMRCellData m_backup;   // Backup solution
  EBAMRCellData m_scratch;  // Scratch data
  EBAMRCellData m_previous; // Previous solution (for each substep)
  EBAMRCellData m_error;    // Error estimator
  EBAMRCellData m_gradient; // Storage for holding gradient
  Vector<EBAMRCellData*> m_extra_storage;

  EBAMRIVData m_scratchIV1; // Scratch data
  EBAMRIVData m_scratchIV2; // Scratch data
  EBAMRIVData m_scratchIV3; // Scratch data
  EBAMRIVData m_scratchIV4; // Scratch data

  EBAMRIFData m_scratchIF1;
  EBAMRIFData m_scratchIF2;
  EBAMRIFData m_scratchIF3;
  EBAMRIFData m_scratchIF4;
};

class strang2::FieldStorage {
public:

  FieldStorage();
  FieldStorage(const int m_stages, const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~FieldStorage();

  void allocate_storage();
  void deallocate_storage();
  
  MFAMRCellData& get_cache(){   return m_cache;   }
  MFAMRCellData& get_backup(){  return m_backup; }
  EBAMRCellData& getElectricFieldCell(){  return m_E_cell;  }
  EBAMRFluxData& getElectricFieldFace(){  return m_E_face;  }
  EBAMRIVData& getElectricFieldEb(){      return m_E_eb;    }
  EBAMRIFData& getElectricFieldDomain(){  return m_E_dom;   }

protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages

  MFAMRCellData m_cache;         // Cached solution
  MFAMRCellData m_backup;        // Scratch storage
  EBAMRCellData m_E_cell;        // Cell-centered E
  EBAMRFluxData m_E_face;        // Face-centered E
  EBAMRIVData m_E_eb;            // EB-centered E
  EBAMRIFData m_E_dom;           // E-field on domain faces
};

class strang2::RtStorage {
public:

  RtStorage();
  RtStorage(const int a_stages, const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~RtStorage();

  void allocate_storage();
  void deallocate_storage();
  
  EBAMRCellData& get_cache(){     return m_cache; }
  EBAMRCellData& get_backup(){    return m_backup; }
  EBAMRIVData& getEbFlux(){     return m_scratchIV; }
  EBAMRIFData& getDomainFlux(){ return m_scratchIF; }

protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages


  EBAMRCellData m_cache;   // Cached solution
  EBAMRCellData m_backup;  // Cached solution
  EBAMRIVData m_scratchIV; // Used for EB-centered isotropic density
  EBAMRIFData m_scratchIF; // Used for domain-centered isotropic density
};

class strang2::SigmaStorage {
public:

  SigmaStorage();
  SigmaStorage(const int a_stages, const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int m_ncomp);
  ~SigmaStorage();
  
  void allocate_storage();
  void deallocate_storage();
  void allocate_extra_storage(const int a_num_extra);
  void deallocate_extra_storage();
  
  EBAMRIVData& get_cache(){   return m_cache; }
  EBAMRIVData& getScratch(){ return m_scratch; }
  EBAMRIVData& get_backup(){  return m_backup; }
  EBAMRIVData& getPrevious(){return m_previous; }
  EBAMRIVData& getError(){   return m_error; }

  Vector<EBAMRIVData*> get_extra_storage(){return m_extra_storage;}

  
protected:

  RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
  phase::which_phase m_phase;    // on this phase
  int m_ncomp;                   // with this many components
  int m_stages;                  // Number of RK stages
  bool m_has_extra;              // For checking for extra storage

  EBAMRIVData m_cache;           // Cached solution
  EBAMRIVData m_backup;          // Backup solution
  EBAMRIVData m_scratch;         // Scratch data
  EBAMRIVData m_previous;        // Previous solution
  EBAMRIVData m_error;           // Error

  Vector<EBAMRIVData*> m_extra_storage;
#include "CD_NamespaceFooter.H"

#endif
