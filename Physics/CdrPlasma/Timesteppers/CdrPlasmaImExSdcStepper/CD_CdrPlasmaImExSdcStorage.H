/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasmaImExSdcStorage.H
   @brief  Declaration of scratch storage inner classes for CdrPlasmaImExSdcStepper
   @author Robert Marskar
*/

#ifndef CD_CDRPLASMAIMEXSDCSTORAGE_H
#define CD_CDRPLASMAIMEXSDCSTORAGE_H

// Our includes
#include <CD_CdrPlasmaImExSdcStepper.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /**
       @brief Per-CDR-species scratch storage for the IMEX-SDC time stepper.
       @details Holds all AMR data arrays needed to advance a single CDR species through the
       SDC integration: the phi^k solution at each of the p+1 quadrature nodes, the explicit
       advection-reaction operator F_AR^k, the implicit diffusion operator F_D^k, the combined
       operator F^k = F_AR^k + F_D^k, and auxiliary arrays for gradient computation, EB/domain
       state extrapolation, and error estimation.
    */
    class CdrPlasmaImExSdcStepper::CdrStorage
    {
    public:
      /**
         @brief Default constructor. Leaves all members uninitialised; call allocateStorage() after construction.
      */
      CdrStorage();

      /**
         @brief Full constructor.
         @param[in] a_amr   AMR mesh used to allocate data.
         @param[in] a_realm Realm on which storage is allocated.
         @param[in] a_phase Phase on which storage is allocated.
         @param[in] m_ncomp Number of components (1 for scalar CDR species).
      */
      CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
                 const std::string             a_realm,
                 const phase::which_phase      a_phase,
                 const int                     m_ncomp);

      /**
         @brief Destructor.
      */
      ~CdrStorage();

      /**
         @brief Allocate all scratch arrays sized for p+1 SDC nodes.
         @param[in] a_p Number of SDC subintervals.
      */
      void
      allocateStorage(const int a_p);

      /**
         @brief Deallocate all scratch arrays.
      */
      void
      deallocateStorage();

      /**
         @brief Return the first cell-centered scratch array (one component).
         @return Reference to m_scratch.
      */
      EBAMRCellData&
      getScratch()
      {
        return m_scratch;
      }

      /**
         @brief Return the second cell-centered scratch array (one component).
         @return Reference to m_scratch2.
      */
      EBAMRCellData&
      getScratch2()
      {
        return m_scratch2;
      }

      /**
         @brief Return the error estimator array (one component, cell-centered).
         @details Stores |phi^{k+1} - phi^k| after each correction sweep.
         @return Reference to m_error.
      */
      EBAMRCellData&
      getError()
      {
        return m_error;
      }

      /**
         @brief Return the cell-centered gradient storage (SpaceDim components).
         @return Reference to m_gradient.
      */
      EBAMRCellData&
      getGradient()
      {
        return m_gradient;
      }

      /**
         @brief Return the "old solution" array used when subcycling is enabled.
         @return Reference to m_old.
      */
      EBAMRCellData&
      getOld()
      {
        return m_old;
      }

      /**
         @brief Return the SpaceDim-component cell-centered scratch array.
         @details Used for intermediate vector quantities such as velocity-weighted states.
         @return Reference to m_scratchD.
      */
      EBAMRCellData&
      getScratchD()
      {
        return m_scratchD;
      }

      /**
         @brief Return the divergence of the advective flux (cell-centered, one component).
         @return Reference to m_divF.
      */
      EBAMRCellData&
      getDivF()
      {
        return m_divF;
      }

      /**
         @brief Return the phi^k solution at all p+1 SDC nodes.
         @return Reference to m_phi. Entry m holds phi^k at node m.
      */
      Vector<EBAMRCellData>&
      getPhi()
      {
        return m_phi;
      }

      /**
         @brief Return the explicit advection-reaction operator F_AR^k at all p+1 nodes.
         @return Reference to m_FAR. Entry m holds F_AR(phi^k_m).
      */
      Vector<EBAMRCellData>&
      getFAR()
      {
        return m_FAR;
      }

      /**
         @brief Return the implicit diffusion operator F_D^k at all p+1 nodes.
         @return Reference to m_FD. Entry m holds F_D(phi^k_m).
      */
      Vector<EBAMRCellData>&
      getFD()
      {
        return m_FD;
      }

      /**
         @brief Return the combined operator F^k = F_AR^k + F_D^k at all p+1 nodes.
         @return Reference to m_F. Entry m holds F(phi^k_m).
      */
      Vector<EBAMRCellData>&
      getF()
      {
        return m_F;
      }

      /**
         @brief Return the first EB-centered scratch array (one component).
         @details Used for diffusion-related EB quantities.
         @return Reference to m_scratchIVs.
      */
      EBAMRIVData&
      getEbScratch1()
      {
        return m_scratchIVs;
      }

      /**
         @brief Return the SpaceDim-component EB-centered scratch array.
         @return Reference to m_scratchIVD.
      */
      EBAMRIVData&
      getEbScratchD()
      {
        return m_scratchIVD;
      }

      /**
         @brief Return the EB-centered species density (extrapolated from the cell center).
         @return Reference to m_scratchIV1.
      */
      EBAMRIVData&
      getEbState()
      {
        return m_scratchIV1;
      }

      /**
         @brief Return the EB-centered normal drift velocity.
         @return Reference to m_scratchIV2.
      */
      EBAMRIVData&
      getEbVelo()
      {
        return m_scratchIV2;
      }

      /**
         @brief Return the EB-centered boundary flux (outward positive).
         @return Reference to m_scratchIV3.
      */
      EBAMRIVData&
      getEbFlux()
      {
        return m_scratchIV3;
      }

      /**
         @brief Return the EB-centered normal density gradient.
         @return Reference to m_scratchIV4.
      */
      EBAMRIVData&
      getEbGrad()
      {
        return m_scratchIV4;
      }

      /**
         @brief Return the domain-face-centered species density (extrapolated from the cell center).
         @return Reference to m_scratchIF1.
      */
      EBAMRIFData&
      getDomainState()
      {
        return m_scratchIF1;
      }

      /**
         @brief Return the domain-face-centered normal drift velocity.
         @return Reference to m_scratchIF2.
      */
      EBAMRIFData&
      getDomainVelo()
      {
        return m_scratchIF2;
      }

      /**
         @brief Return the domain-face-centered boundary flux (inward positive).
         @return Reference to m_scratchIF3.
      */
      EBAMRIFData&
      getDomainFlux()
      {
        return m_scratchIF3;
      }

      /**
         @brief Return the domain-face-centered normal density gradient.
         @return Reference to m_scratchIF4.
      */
      EBAMRIFData&
      getDomainGrad()
      {
        return m_scratchIF4;
      }

    protected:
      /**
         @brief Realm on which storage is allocated.
      */
      std::string m_realm;

      /**
         @brief AMR mesh used for allocation.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Phase on which storage is allocated.
      */
      phase::which_phase m_phase;

      /**
         @brief Number of scalar components (1 for CDR species).
      */
      int m_ncomp;

      /**
         @brief Number of SDC subintervals (p+1 nodes).
      */
      int m_p;

      /**
         @brief If true, extra subcycling storage is allocated.
      */
      bool m_subcycle;

      /**
         @brief First general-purpose cell-centered scratch array.
      */
      EBAMRCellData m_scratch;

      /**
         @brief Second general-purpose cell-centered scratch array.
      */
      EBAMRCellData m_scratch2;

      /**
         @brief Cell-centered error estimator: |phi^{k+1} - phi^k|.
      */
      EBAMRCellData m_error;

      /**
         @brief Cell-centered gradient of the species density (SpaceDim components).
      */
      EBAMRCellData m_gradient;

      /**
         @brief Cell-centered solution at the start of the step. Used for subcycling.
      */
      EBAMRCellData m_old;

      /**
         @brief SpaceDim-component cell-centered scratch array.
      */
      EBAMRCellData m_scratchD;

      /**
         @brief Cell-centered divergence of the advective flux.
      */
      EBAMRCellData m_divF;

      /**
         @brief phi^k at each of the p+1 SDC nodes. Index m ∈ [0, p].
      */
      Vector<EBAMRCellData> m_phi;

      /**
         @brief Explicit advection-reaction operator F_AR(phi^k) at each node. Index m ∈ [0, p].
      */
      Vector<EBAMRCellData> m_FAR;

      /**
         @brief Implicit diffusion operator F_D(phi^k) at each node. Index m ∈ [0, p].
      */
      Vector<EBAMRCellData> m_FD;

      /**
         @brief Combined operator F = F_AR + F_D at each node. Index m ∈ [0, p].
      */
      Vector<EBAMRCellData> m_F;

      /**
         @brief EB-centered scalar scratch array (one component, e.g. for diffusion states).
      */
      EBAMRIVData m_scratchIVs;

      /**
         @brief EB-centered SpaceDim-component scratch array.
      */
      EBAMRIVData m_scratchIVD;

      /**
         @brief EB-centered species density extrapolated from the cell center.
      */
      EBAMRIVData m_scratchIV1;

      /**
         @brief EB-centered normal drift velocity component.
      */
      EBAMRIVData m_scratchIV2;

      /**
         @brief EB-centered boundary flux.
      */
      EBAMRIVData m_scratchIV3;

      /**
         @brief EB-centered normal density gradient.
      */
      EBAMRIVData m_scratchIV4;

      /**
         @brief Domain-face-centered species density extrapolated from the cell center.
      */
      EBAMRIFData m_scratchIF1;

      /**
         @brief Domain-face-centered normal drift velocity component.
      */
      EBAMRIFData m_scratchIF2;

      /**
         @brief Domain-face-centered boundary flux.
      */
      EBAMRIFData m_scratchIF3;

      /**
         @brief Domain-face-centered normal density gradient.
      */
      EBAMRIFData m_scratchIF4;
    };

    /**
       @brief Field solver scratch storage for the IMEX-SDC time stepper.
       @details Holds the electrostatic potential at the previous step, and the electric field
       at all locations needed by the CDR and surface charge solvers: cell centers, face centers,
       EB centroids, and domain faces.
    */
    class CdrPlasmaImExSdcStepper::FieldStorage
    {
    public:
      /**
         @brief Default constructor. Leaves all members uninitialised; call allocateStorage() after construction.
      */
      FieldStorage();

      /**
         @brief Full constructor.
         @param[in] a_amr    AMR mesh used to allocate data.
         @param[in] a_realm  Realm on which storage is allocated.
         @param[in] a_phase  Phase on which storage is allocated.
         @param[in] a_ncomp  Number of components (SpaceDim for the electric field).
      */
      FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
                   const std::string             a_realm,
                   const phase::which_phase      a_phase,
                   const int                     a_ncomp);

      /**
         @brief Destructor.
      */
      ~FieldStorage();

      /**
         @brief Allocate all field storage arrays for p+1 SDC nodes.
         @param[in] a_p Number of SDC subintervals.
      */
      void
      allocateStorage(const int a_p);

      /**
         @brief Deallocate all field storage arrays.
      */
      void
      deallocateStorage();

      /**
         @brief Return the multifluid electrostatic potential at the previous time step.
         @return Reference to m_previous.
      */
      MFAMRCellData&
      getPrevious()
      {
        return m_previous;
      }

      /**
         @brief Return the cell-centered electric field vector (SpaceDim components).
         @return Reference to m_E_cell.
      */
      EBAMRCellData&
      getElectricFieldCell()
      {
        return m_E_cell;
      }

      /**
         @brief Return the face-centered normal electric field component.
         @return Reference to m_E_face.
      */
      EBAMRFluxData&
      getElectricFieldFace()
      {
        return m_E_face;
      }

      /**
         @brief Return the EB-centroid-centered electric field vector.
         @return Reference to m_E_eb.
      */
      EBAMRIVData&
      getElectricFieldEb()
      {
        return m_E_eb;
      }

      /**
         @brief Return the domain-face-centered electric field vector.
         @return Reference to m_E_dom.
      */
      EBAMRIFData&
      getElectricFieldDomain()
      {
        return m_E_dom;
      }

    protected:
      /**
         @brief Realm on which storage is allocated.
      */
      std::string m_realm;

      /**
         @brief AMR mesh used for allocation.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Phase on which storage is allocated.
      */
      phase::which_phase m_phase;

      /**
         @brief Number of field components (SpaceDim).
      */
      int m_ncomp;

      /**
         @brief Number of SDC subintervals.
      */
      int m_p;

      /**
         @brief Multifluid cell-centered electrostatic potential at the previous time step.
      */
      MFAMRCellData m_previous;

      /**
         @brief Cell-centered electric field vector (SpaceDim components).
      */
      EBAMRCellData m_E_cell;

      /**
         @brief Face-centered normal component of the electric field.
      */
      EBAMRFluxData m_E_face;

      /**
         @brief EB-centroid-centered electric field vector (SpaceDim components).
      */
      EBAMRIVData m_E_eb;

      /**
         @brief Domain-face-centered electric field vector (SpaceDim components).
      */
      EBAMRIFData m_E_dom;
    };

    /**
       @brief Per-RTE-solver scratch storage for the IMEX-SDC time stepper.
       @details Holds the photon density at the previous time step plus EB and domain face
       flux arrays used to impose boundary conditions during the RTE advance.
    */
    class CdrPlasmaImExSdcStepper::RtStorage
    {
    public:
      /**
         @brief Default constructor. Leaves all members uninitialised; call allocateStorage() after construction.
      */
      RtStorage();

      /**
         @brief Full constructor.
         @param[in] a_amr    AMR mesh used to allocate data.
         @param[in] a_realm  Realm on which storage is allocated.
         @param[in] a_phase  Phase on which storage is allocated.
         @param[in] a_ncomp  Number of components (1 for scalar RTE species).
      */
      RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
                const std::string             a_realm,
                const phase::which_phase      a_phase,
                const int                     a_ncomp);

      /**
         @brief Destructor.
      */
      ~RtStorage();

      /**
         @brief Allocate RTE scratch storage for p+1 SDC nodes.
         @param[in] a_p Number of SDC subintervals.
      */
      void
      allocateStorage(const int a_p);

      /**
         @brief Deallocate all RTE scratch storage.
      */
      void
      deallocateStorage();

      /**
         @brief Return the cell-centered photon density at the previous time step.
         @return Reference to m_previous.
      */
      EBAMRCellData&
      getPrevious()
      {
        return m_previous;
      }

      /**
         @brief Return the EB-centered isotropic photon flux scratch array.
         @return Reference to m_scratchIV.
      */
      EBAMRIVData&
      getEbFlux()
      {
        return m_scratchIV;
      }

      /**
         @brief Return the domain-face-centered isotropic photon flux scratch array.
         @return Reference to m_scratchIF.
      */
      EBAMRIFData&
      getDomainFlux()
      {
        return m_scratchIF;
      }

    protected:
      /**
         @brief Realm on which storage is allocated.
      */
      std::string m_realm;

      /**
         @brief AMR mesh used for allocation.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Phase on which storage is allocated.
      */
      phase::which_phase m_phase;

      /**
         @brief Number of scalar components (1 for RTE species).
      */
      int m_ncomp;

      /**
         @brief Number of SDC subintervals.
      */
      int m_p;

      /**
         @brief Cell-centered photon density at the previous time step (used to restore state on rejection).
      */
      EBAMRCellData m_previous;

      /**
         @brief EB-centered scratch array for the isotropic photon density boundary condition.
      */
      EBAMRIVData m_scratchIV;

      /**
         @brief Domain-face-centered scratch array for the isotropic photon density boundary condition.
      */
      EBAMRIFData m_scratchIF;
    };

    /**
       @brief Surface charge solver scratch storage for the IMEX-SDC time stepper.
       @details Holds sigma^k at each SDC node plus the new and lagged flux operator values
       needed for the SDC correction sweeps, and an error estimator.
    */
    class CdrPlasmaImExSdcStepper::SigmaStorage
    {
    public:
      /**
         @brief Default constructor. Leaves all members uninitialised; call allocateStorage() after construction.
      */
      SigmaStorage();

      /**
         @brief Full constructor.
         @param[in] a_amr    AMR mesh used to allocate data.
         @param[in] a_realm  Realm on which storage is allocated.
         @param[in] a_phase  Phase on which storage is allocated.
         @param[in] a_ncomp  Number of components (1 for scalar surface charge).
      */
      SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
                   const std::string             a_realm,
                   const phase::which_phase      a_phase,
                   const int                     a_ncomp);

      /**
         @brief Destructor.
      */
      ~SigmaStorage();

      /**
         @brief Allocate surface charge scratch storage for p+1 SDC nodes.
         @param[in] a_p Number of SDC subintervals.
      */
      void
      allocateStorage(const int a_p);

      /**
         @brief Deallocate all surface charge scratch storage.
      */
      void
      deallocateStorage();

      /**
         @brief Return the general-purpose EB-centered scratch array.
         @return Reference to m_scratch.
      */
      EBAMRIVData&
      getScratch()
      {
        return m_scratch;
      }

      /**
         @brief Return the surface charge error estimator: |sigma^{k+1} - sigma^k|.
         @return Reference to m_error.
      */
      EBAMRIVData&
      getError()
      {
        return m_error;
      }

      /**
         @brief Return sigma^k at all p+1 SDC nodes.
         @return Reference to m_sigma. Entry m holds sigma^k at node m.
      */
      Vector<EBAMRIVData>&
      getSigmaSolver()
      {
        return m_sigma;
      }

      /**
         @brief Return the new (current-sweep) flux operator slopes at all p+1 nodes.
         @return Reference to m_Fnew. Entry m holds F_{sigma,new}^k at node m.
      */
      Vector<EBAMRIVData>&
      getFnew()
      {
        return m_Fnew;
      }

      /**
         @brief Return the lagged (previous-sweep) flux operator slopes at all p+1 nodes.
         @return Reference to m_Fold. Entry m holds F_{sigma,old}^{k-1} at node m.
      */
      Vector<EBAMRIVData>&
      getFold()
      {
        return m_Fold;
      }

    protected:
      /**
         @brief Realm on which storage is allocated.
      */
      std::string m_realm;

      /**
         @brief AMR mesh used for allocation.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Phase on which storage is allocated.
      */
      phase::which_phase m_phase;

      /**
         @brief Number of scalar components (1 for surface charge).
      */
      int m_ncomp;

      /**
         @brief Number of SDC subintervals.
      */
      int m_p;

      /**
         @brief General-purpose EB-centered scratch array.
      */
      EBAMRIVData m_scratch;

      /**
         @brief Surface charge error estimator: |sigma^{k+1} - sigma^k|.
      */
      EBAMRIVData m_error;

      /**
         @brief sigma^k at each of the p+1 SDC nodes. Index m ∈ [0, p].
      */
      Vector<EBAMRIVData> m_sigma;

      /**
         @brief New (current-sweep) flux operator slopes at each node. Index m ∈ [0, p].
      */
      Vector<EBAMRIVData> m_Fnew;

      /**
         @brief Lagged (previous-sweep) flux operator slopes at each node. Index m ∈ [0, p].
      */
      Vector<EBAMRIVData> m_Fold;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
