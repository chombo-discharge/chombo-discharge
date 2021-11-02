/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStorage.H
  @brief  Declaration of scratch storage for CdrPlasmaGodunovStepper
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaGodunovStorage_H
#define CD_CdrPlasmaGodunovStorage_H

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
  
    class CdrPlasmaGodunovStepper::CdrStorage {
    public:

      CdrStorage();
      CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
		 const std::string              a_realm,
		 const phase::which_phase       a_phase,
		 const int m_ncomp);
    
      ~CdrStorage();

      void allocateStorage();
      void deallocateStorage();

      EBAMRCellData& getScratch()  { return m_scratch;    }
      EBAMRCellData& getScratch2() { return m_scratch2;   }
      EBAMRCellData& getScratch3() { return m_scratch3;   }
      EBAMRCellData& getExtrap()   { return m_cellExtr;   }
      EBAMRCellData& getGradient() { return m_gradient;   }
      EBAMRCellData& getDivD()     { return m_divDgradPhi;}      

      EBAMRIVData& getEbScratch1(){ return m_scratchIVs; }
      EBAMRIVData& getEbScratchD(){ return m_scratchIVD; }
      EBAMRIVData& getEbState(){    return m_scratchIV1; }
      EBAMRIVData& getEbVelo(){     return m_scratchIV2; }
      EBAMRIVData& getEbFlux(){     return m_scratchIV3; }
      EBAMRIVData& getEbGrad(){     return m_scratchIV4; }

      EBAMRIFData& getDomainState(){  return m_scratchIF1; }
      EBAMRIFData& getDomainVelo(){   return m_scratchIF2; }
      EBAMRIFData& getDomainFlux(){   return m_scratchIF3; }
      EBAMRIFData& getDomainGrad(){   return m_scratchIF4; }
  
    protected:
      std::string             m_realm;
      RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
      phase::which_phase m_phase;    // on this phase
      int m_ncomp;                   // with this many components

      EBAMRCellData m_scratch;     // Scratch data
      EBAMRCellData m_scratch2;    // Scratch data
      EBAMRCellData m_scratch3;    // Scratch data
      EBAMRCellData m_cellExtr;    // Scratch data
      EBAMRCellData m_gradient;    // Storage for holding gradient
      EBAMRCellData m_divDgradPhi; // Storage for holding div(D*phi) when using semi-ipmlicit formulation.

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

    class CdrPlasmaGodunovStepper::FieldStorage {
    public:

      FieldStorage();
      FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
		   const std::string              a_realm,
		   const phase::which_phase       a_phase,
		   const int m_ncomp);
      ~FieldStorage();

      void allocateStorage();
      void deallocateStorage();

      EBAMRCellData& getElectricFieldCell()      {return m_E_cell;   }
      EBAMRIVData& getElectricFieldEb()          {return m_E_eb;     }
      EBAMRIFData& getElectricFieldDomain()      {return m_E_dom;    }

    protected:
      std::string m_realm;
      RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
      phase::which_phase m_phase;    // on this phase
      int m_ncomp;                   // with this many components

      EBAMRCellData m_E_cell;        // Cell-centered E
      EBAMRIVData m_E_eb;            // EB-centered E
      EBAMRIFData m_E_dom;           // E-field on domain faces
    };

    class CdrPlasmaGodunovStepper::RtStorage {
    public:

      RtStorage();
      RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
		const std::string              a_realm,
		const phase::which_phase       a_phase,
		const int                      a_ncomp);
      ~RtStorage();

      void allocateStorage();
      void deallocateStorage();

      EBAMRIVData& getEbFlux(){     return m_scratchIV; }
      EBAMRIFData& getDomainFlux(){ return m_scratchIF; }

    protected:

      std::string m_realm;
      RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
      phase::which_phase m_phase;    // on this phase
      int m_ncomp;                   // with this many components

      EBAMRIVData m_scratchIV;
      EBAMRIFData m_scratchIF;
    };

    class CdrPlasmaGodunovStepper::SigmaStorage {
    public:

      SigmaStorage();
      SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
		   const std::string              a_realm,
		   const phase::which_phase       a_phase,
		   const int                      a_ncomp);
      ~SigmaStorage();
  
      void allocateStorage();
      void deallocateStorage();

      EBAMRIVData& getScratch() {return m_scratch;}
  
    protected:

      std::string m_realm;
      RefCountedPtr<AmrMesh> m_amr; // Needed for allocation. Will allocate
      phase::which_phase m_phase;    // on this phase
      int m_ncomp;                   // with this many components

      EBAMRIVData m_scratch;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
