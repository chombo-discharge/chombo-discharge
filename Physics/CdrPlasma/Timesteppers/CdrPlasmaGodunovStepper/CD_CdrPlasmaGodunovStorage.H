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

    /*!
      @brief Utility class for memory handling in CdrPlasmaStepper.
      @details This class exists because we often need some extra storage when advancing the equations. 
    */
    class CdrPlasmaGodunovStepper::CdrStorage {
    public:

      /*!
	@brief Disallowed constructor. Only strong construction allowed.
      */
      CdrStorage() = delete;

      /*!
	@brief Full constructor. After this, user may call allocateStorage and deallocateStorage
	@param[in] a_amr   Reference to AmrMesh
	@param[in] a_realm Realm where storage lives
	@param[in] a_phase Phase where storage lives (only single-phased data in here)
      */
      CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
		 const std::string             a_realm,
		 const phase::which_phase      a_phase);

      /*!
	@brief Destructor
      */
      virtual ~CdrStorage();

      /*!
	@brief Allocation function. This will allocate all transient storage in this class
      */
      virtual void allocateStorage();

      /*!
	@brief Allocation function. This will delete all transient storage in this class
      */      
      virtual void deallocateStorage();

      /*!
	@brief Get scratch storage
      */
      EBAMRCellData& getScratch() const {
	return m_scratch;
      }

      /*!
	@brief Get the other scratch storage. 
      */      
      EBAMRCellData& getScratch2() const {
	return m_scratch2;
      }

      /*!
	@brief Get the extrapolation storage. This is usually used for holding phi + 0.5*dt*S
      */
      EBAMRCellData& getExtrap() const {
	return m_cellExtr;
      }

      /*!
	@brief Get storage for holding gradient. 
      */      
      EBAMRCellData& getGradient() const {
	return m_gradient;
      }

      /*!
	@brief Get scratch storage that holds the density on the EB
      */
      EBAMRIVData& getEbState() const {
	return m_scratchIV1;
      }

      /*!
	@brief Get the scratch storage holding the normal component of the velocity on the EB
      */
      EBAMRIVData& getEbVelo() const {
	return m_scratchIV2;
      }

      /*!
	@brief Get the scratch storage holding the extrapolated flux on the EB
      */      
      EBAMRIVData& getEbFlux() const {
	return m_scratchIV3;
      }

      /*!
	@brief Get the scratch storage holding the extrapolated gradient on the EB
      */            
      EBAMRIVData& getEbGrad() const {
	return m_scratchIV4;
      }

      /*!
	@brief Get the scratch storage holding the extrapolated density on domain faces
      */
      EBAMRIFData& getDomainState() const {
	return m_scratchIF1;
      }

      /*!
	@brief Get the scratch storage holding the extrapolated velocity (normal component only) on domain faces
      */      
      EBAMRIFData& getDomainVelo() const {
	return m_scratchIF2;
      }

      /*!
	@brief Get the scratch storage holding the extrapolated flux on domain faces
      */            
      EBAMRIFData& getDomainFlux() const {
	return m_scratchIF3;
      }

      /*!
	@brief Get the scratch storage holding the extrapolated gradient (normal component only) on domain faces
      */                  
      EBAMRIFData& getDomainGrad() const {
	return m_scratchIF4;
      }
  
    protected:

      /*!
	@brief Reference to AmrMesh
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Realm where data lives
      */
      std::string m_realm;

      /*!
	@brief Phase where data lives
      */      
      phase::which_phase m_phase;

      /*!
	@brief Cell-centered scratch data
      */
      mutable EBAMRCellData m_scratch;     

      /*!
	@brief Cell-centered scratch data
      */      
      mutable EBAMRCellData m_scratch2;

      /*!
	@brief Cell-centered scratch data. Used for holding phi + 0.5*dt*S
      */
      mutable EBAMRCellData m_cellExtr;

      /*!
	@brief Cell-centered scratch data for holding the gradient.
      */      
      mutable EBAMRCellData m_gradient;

      /*!
	@brief EB-centered scratch data
      */
      mutable EBAMRIVData m_scratchIV1;

      /*!
	@brief EB-centered scratch data
      */      
      mutable EBAMRIVData m_scratchIV2;

      /*!
	@brief EB-centered scratch data
      */      
      mutable EBAMRIVData m_scratchIV3;

      /*!
	@brief EB-centered scratch data
      */      
      mutable EBAMRIVData m_scratchIV4;

      /*!
	@brief Domain-centered scratch data
      */            
      mutable EBAMRIFData m_scratchIF1;
      /*!
	@brief Domain-centered scratch data
      */                  
      mutable EBAMRIFData m_scratchIF2;

      /*!
	@brief Domain-centered scratch data
      */                  
      mutable EBAMRIFData m_scratchIF3;

      /*!
	@brief Domain-centered scratch data
      */                  
      mutable EBAMRIFData m_scratchIF4;
    };


    /*!
      @brief Utility class for memory handling in CdrPlasmaStepper.
      @details This class exists because we often need some extra storage when advancing the equations. 
    */    
    class CdrPlasmaGodunovStepper::FieldStorage {
    public:

      /*!
	@brief Disallowed weak constructor
      */
      FieldStorage() = delete;

      /*!
	@brief Full constructor. Can call allocation functions after this. 
	@param[in] a_amr   Reference to AmrMesh
	@param[in] a_realm Realm where storage lives
	@param[in] a_phase Phase where storage lives (only single-phased data in here)
      */      
      FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
		   const std::string             a_realm,
		   const phase::which_phase      a_phase);

      /*!
	@brief Destructor
      */
      virtual ~FieldStorage();

      /*!
	@brief Allocation function -- allocates all the necessary memory. 
      */
      virtual void allocateStorage();

      /*!
	@brief Deallocation function -- deallocates all the necessary memory. 
      */      
      virtual void deallocateStorage();

      /*!
	@brief Get scratch storage that holds the electric field at the cell center
      */
      EBAMRCellData& getElectricFieldCell() const {
	return m_electricFieldCell;
      }

      /*!
	@brief Get scratch storage that holds the electric field at the EB
      */      
      EBAMRIVData& getElectricFieldEB() const {
	return m_electricFieldEB;
      }

      /*!
	@brief Get scratch storage that holds the electric field at the domainf aces
      */      
      EBAMRIFData& getElectricFieldDomain() const {
	return m_electricFieldDomain;
      }

    protected:

      /*!
	@brief Reference to AmrMesh
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Realm where data lives
      */
      std::string m_realm;

      /*!
	@brief Phase where data lives
      */      
      phase::which_phase m_phase;      

      /*!
	@brief Scratch storage for holding electric field at cell center.
      */
      mutable EBAMRCellData m_electricFieldCell;

      /*!
	@brief Scratch storage for holding electric field at the EB.
      */      
      mutable EBAMRIVData m_electricFieldEB;

      /*!
	@brief Scratch storage for holding electric field at domain faces. 
      */            
      mutable EBAMRIFData m_electricFieldDomain;
    };

    /*!
      @brief Utility class for memory handling in CdrPlasmaStepper.
      @details This class exists because we often need some extra storage when advancing the equations. 
    */        
    class CdrPlasmaGodunovStepper::RtStorage {
    public:

      /*!
	@brief Disallowed weak constructor
      */
      RtStorage() = delete;

      /*!
	@brief Full constructor. Can call allocation functions after this. 
	@param[in] a_amr   Reference to AmrMesh
	@param[in] a_realm Realm where storage lives
	@param[in] a_phase Phase where storage lives (only single-phased data in here)
      */      
      RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
		const std::string             a_realm,
		const phase::which_phase      a_phase);

      /*!
	@brief Destructor
      */
      virtual ~RtStorage();

      /*!
	@brief Allocate function -- allocates the flux data holders
      */
      virtual void allocateStorage();

      /*!
	@brief Allocate function -- deallocates the flux data holders
      */      
      virtual void deallocateStorage();

      /*!
	@brief Get scratch storage that holds the EB flux
      */
      EBAMRIVData& getEbFlux() const {
	return m_scratchIV;
      }

      /*!
	@brief Get scratch storage that holds the domain
      */      
      EBAMRIFData& getDomainFlux() const {
	return m_scratchIF;
      }

    protected:

      /*!
	@brief Reference to AmrMesh
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Realm where data lives
      */
      std::string m_realm;

      /*!
	@brief Phase where data lives
      */      
      phase::which_phase m_phase;            

      /*!
	@brief Scratch storage on EB
      */
      mutable EBAMRIVData m_scratchIV;

      /*!
	@brief Scratch storage on domain faces
      */      
      mutable EBAMRIFData m_scratchIF;
    };

    /*!
      @brief Utility class for memory handling in CdrPlasmaStepper.
      @details This class exists because we often need some extra storage when advancing the equations. 
    */            
    class CdrPlasmaGodunovStepper::SigmaStorage {
    public:

      /*!
	@brief Disallowed weak constructor. Use full construction.
      */
      SigmaStorage() = delete;

      /*!
	@brief Full constructor. Can call allocation functions after this. 
	@param[in] a_amr   Reference to AmrMesh
	@param[in] a_realm Realm where storage lives
	@param[in] a_phase Phase where storage lives (only single-phased data in here)
      */            
      SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
		   const std::string             a_realm,
		   const phase::which_phase      a_phase);

      /*!
	@brief Destructor
      */
      virtual ~SigmaStorage();

      /*!
	@brief Allocation function -- allocates the storage
      */
      virtual void allocateStorage();

      /*!
	@brief Allocation function -- deallocates the storage
      */      
      virtual void deallocateStorage();

      /*!
	@brief Get scratch data
      */
      EBAMRIVData& getScratch() const {
	return m_scratch;
      }
  
    protected:

      /*!
	@brief Reference to AmrMesh
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Realm where data lives
      */
      std::string m_realm;

      /*!
	@brief Phase where data lives
      */      
      phase::which_phase m_phase;            

      /*!
	@brief Scratch data
      */
      mutable EBAMRIVData m_scratch;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
