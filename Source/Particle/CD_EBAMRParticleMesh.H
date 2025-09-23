/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRParticleMesh.H
  @brief  Declaration of a class for handling particle-mesh operations with AMR.
  @author Robert Marskar
*/

#ifndef CD_EBAMRParticleMesh_H
#define CD_EBAMRParticleMesh_H

// Chombo includes
#include <RefCountedPtr.H>
#include <Copier.H>
#include <EBLevelGrid.H>

// Our includes
#include <CD_DepositionType.H>
#include <CD_EBAMRData.H>
#include <CD_CoarseFineDeposition.H>
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBParticleMesh.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for handling particle-mesh operations with AMR. 
  @details This class contains templated functions for depositing particles with AMR. When depositing, the user must specify the deposition type (e.g., CIC) as
  well as how to handle coarse-fine boundary conditions and cut-cell deposition. Because of the many degrees of freedom in depositing particles with AMR (and EB), 
  the user must expect to preprocess his particle data before solving. There are several ways deposition works through this class:

  When specifying CoarseFineDeposition::Interp, the deposition is done by first depositing the data as usual on both levels. On the fine level, there may be mass that
  hangs over the refinement boundaries; that mass is added to the coarse level. On the coarse level some of the particles will deposit mass into cells that are covered
  by a finer level; that mass is interpolated with piece-wise constant interpolation and then added to the fine level. 

  When specifying CoarseFineDeposition::Halo, the handling of the coarse-side deposition is different. The mass from the fine-level particles that hangs over the 
  refinement boundaries is put on the coarse mesh. However, the coarse-level particles are deposited in two separate steps so that the coarse-level particles
  have the same physical widths on the coarse and fine levels. We first deposit ALL the coarse-level particles on the coarse level. Some of these particles (halo particles)
  will have clouds that overlap with the fine level. We take the set of these particles and deposit them directly onto the fine level (via buffers) with 2x or 4x the
  original particle width (depending on the refinement ratios between the levels). Essentially, this differs CoarseFineDeposition::Interp in the sense that the original 
  particle width is maintained on the fine level.

  When specifying CoarseFineDeposition::HaloNGP, the particles whose deposition clouds hang over the refinement boundary are deposited on the fine level, and the ghosted
  mass is then added to the coarse side. Particles that live on the coarse side of a refinement boundary are deposited using an NGP scheme.

  When specifying CoarseFineDeposition::Transition, the particles immediately on the coarse-side of the refinement boundary deposit with the fine-grid particle width. 
*/
class EBAMRParticleMesh
{
public:
  /*!
    @brief Default constructor. Leaves object in undefined state and thus requires the user to call define.
  */
  EBAMRParticleMesh();

  /*!
    @brief Full constructor. Calls the define function
    @param[in] a_eblgs            Grids
    @param[in] a_refRat           Refinement factory between levels
    @param[in] a_dx               Grid resolutions.
    @πaram[in] a_probLo           Lower-left corner of physical domain.
    @param[in] a_ghost            Ghost cells in data holders.
    @param[in] a_finestLevel      Finest grid level
  */
  EBAMRParticleMesh(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                    const Vector<int>&                        a_refRat,
                    const Vector<Real>&                       a_dx,
                    const RealVect&                           a_probLo,
                    const int&                                a_ghost,
                    const int                                 a_finestLevel);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBAMRParticleMesh();

  /*!
    @brief Define function
    @param[in] a_eblgs            Grids
    @param[in] a_refRat           Refinement factory between levels
    @param[in] a_dx               Grid resolutions.
    @πaram[in] a_probLo           Lower-left corner of physical domain.
    @param[in] a_ghost            Ghost cells in data holders.
    @param[in] a_finestLevel      Finest grid level
  */
  void
  define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
         const Vector<int>&                        a_refRat,
         const Vector<Real>&                       a_dx,
         const RealVect&                           a_probLo,
         const int&                                a_ghost,
         const int                                 a_finestLevel);

  /*!
    @brief Class for deposition of particles of a type P to the mesh. This is the main function that users should call. 
    @details P is the particle type, Ret is the returned value of the particle member function, and must be a Real or a
    RealVect. MemFunc is a pointer to a member function of P. E.g., deposit<P, const RealVect&, &P::position>.
    @param[out] a_meshData             Mesh data. Must have exactly one compnent. 
    @πaram[in]  a_particles            Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType       Specification of deposition kernel (e.g., CIC)
    @param[in]  a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]  a_forceIrregNGP        Force NGP deposition in irregular cells or not. 
    @note When using certain types of coarseFineDeposition methods -- the user may need to fill the mask particles before calling this routine. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)() const>
  void
  deposit(EBAMRCellData&              a_meshData,
          const ParticleContainer<P>& a_particles,
          const DepositionType        a_depositionType,
          const CoarseFineDeposition  a_coarseFineDeposition,
          const bool                  a_forceIrregNGP = false);

  /*!
    @brief Interpolate a scalar field onto the particle position. 
    @details This is just like regular particle-mesh interpolation. The input field should have exactly one component and the
    the field will be interpolated onto the template parameter's input field. The template parameter 'particleScalarField' should be a pointer
    to a member function which will set the particle field. E.g. the function must have a signature Real& P::particleScalarField(). A valid
    expression is e.g.
    
    interpolate<P, &P::mass> (...)
	
    To interpolate onto a different field, replace ::mass by the other scalar field.

    @param[inout] a_particles       Particles to be interpolated. 
    @param[in]    a_meshScalarField Scalar field on the mesh 
    @param[in]    a_interpType      Interpolation type. 
    @param[in]    a_forceIrregNGP   Force NGP interpolation in cut-cells. 
  */
  template <class P, Real& (P::*particleScalarField)()>
  void
  interpolate(ParticleContainer<P>& a_particles,
              const EBAMRCellData&  a_meshVectorField,
              const DepositionType  a_interpType,
              const bool            a_forceIrregNGP = false) const;

  /*!
    @brief Interpolate a vector field onto the particle position. 
    @details This is just like regular particle-mesh interpolation. The input field should have exactly SpaceDim components and the
    the field will be interpolated onto the template parameter's input field. The template parameter 'particleVectorField' should be a pointer
    to a member function which will set the particle field. E.g. the function must have a signature RealVect& P::particleVectorField(). A valid
    expression is e.g.
    
    interpolate<P, &P::velocity> (...)
	
    To interpolate onto a different field, replace ::velocity by another vector field. 

    @param[inout] a_particles            Particles to be interpolated. 
    @param[in]    a_meshScalarField      Scalar field on the mesh 
    @param[in]    a_interpType           Interpolation type. 
    @param[in]    a_coarseFineDeposition Specification of handling of coarse-fine boundaries.
    @param[in]    a_forceIrregNGP        Force NGP interpolation in cut-cells. 
  */
  template <class P, RealVect& (P::*particleVectorField)()>
  void
  interpolate(ParticleContainer<P>& a_particles,
              const EBAMRCellData&  a_meshVectorField,
              const DepositionType  a_interpType,
              const bool            a_forceIrregNGP = false) const;

  /*!
    @brief Get buffers for handling deposition over refinement boundaries.
    @return The buffer for handling mass transfer between level lvl and level lvl+1 lives on lvl+1. 
  */
  Vector<RefCountedPtr<EBCoarseFineParticleMesh>>&
  getEBCoarseFineParticleMesh() const;

  /*!
    @brief Get EBParticleMesh deposition/interpolation object for specified grid patch.
    @details This routine is used for when users want to interpolate on a per-patch basis. This is useful when the user has written his interpolation
    functions on a per-patch basis. I fail to see any other reason where this routine is useful. 
    @param[in] a_lvl Grid level
    @param[in] a_dit Grid index
  */
  const EBParticleMesh&
  getEBParticleMesh(const int a_lvl, const DataIndex& a_dit) const;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Verbose or not
  */
  bool m_verbose;

  /*!
    @brief Lower-left corner of physical domain
  */
  RealVect m_probLo;

  /*!
    @brief Number of ghost cells
  */
  int m_ghost;

  /*!
    @brief Finest AMR level
  */
  int m_finestLevel;

  /*!
    @brief Grids on each level
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgs;

  /*!
    @brief Refinement ratios between levels
  */
  Vector<int> m_refRat;

  /*!
    @brief Grid resolutions.
  */
  Vector<Real> m_dx;

  /*!
    @brief Regular particle-mesh object on each grid level
  */
  Vector<RefCountedPtr<LayoutData<EBParticleMesh>>> m_ebParticleMesh;

  /*!
    @brief Special particle-mesh objects for depositing on the coarsened fine grid. 
    @details These are deposition/interpolation objects defined on the refined coarse grid. They are used when we need
    to deposit coarse-grid particles on the fine grid. 
  */
  Vector<RefCountedPtr<LayoutData<EBParticleMesh>>> m_ebParticleMeshFiCo;

  /*!
    @brief Valid grid cells. 
  */
  Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> m_validCells;

  /*!
    @brief Outer halo masks with various widths.
    @details These masks are true on the coarse side on the refinement boundary.
    @note We have purposefully undefined the mask in patches where the mask is always false. This means
    that one should never call the copy function with this mask.
  */
  std::map<int, Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>> m_outerHaloMasks;

  /*!
    @brief Transition masks with various widths.
    @details These masks are true in the fine-grid ghost cells on the coarse-side of the refinement boundary.
    @note We have purposefully undefined the mask in patches where the mask is always false. This means
    that one should never call the copy function with this mask.    
  */
  std::map<int, Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>> m_transitionMasks;

  /*!
    @brief Buffers for handling arithmetic for mass moving from coarse to fine level and vice versa. 
    @note The buffer for handling mass transfer between level lvl and level lvl+1 lives on lvl+1. 
  */
  mutable Vector<RefCountedPtr<EBCoarseFineParticleMesh>> m_coarseFinePM;

  /*!
    @brief Copier for moving data from valid+ghost to valid on each AMR level.
    @note After depositing on a level, this is used in order to add data from the ghost cells to the valid region. 
  */
  mutable Vector<Copier> m_levelCopiers;

  /*!
    @brief Define level copiers
    @note This defines the Copiers we need to move data from valid+ghost to valid on each level.
  */
  void
  defineLevelMotion();

  /*!
    @brief Define coarse-fine data motion operators
    @details This defines the necessary constructors for 
  */
  void
  defineCoarseFineMotion();

  /*!
    @brief Define EBParticleMesh objects
  */
  void
  defineEBParticleMesh();

  /*!
    @brief Define the outer halo masks
  */
  void
  defineOuterHaloMasks();

  /*!
    @brief Define the coarse-fine transition masks
  */
  void
  defineTransitionMasks();

  /*!
    @brief Get transition mask width
    @details Used when depositing with a transition zone.
    @param[in] a_depositionType Deposition type
    @param[in] a_refRat Refinement ratio between levels
  */
  int
  getTransitionMaskWidth(const DepositionType a_depositionType, const int a_refRat) const;

  /*!
    @brief Deposit particles on the mesh, interpolating the coarse-grid invalid mass to the fine grid.     
    @param[out] a_meshData       Mesh data. 
    @param[in]  a_particles      Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType Specification of deposition kernel (e.g., CIC)
    @param[in]  a_forceIrregNGP  Force NGP deposition in irregular cells or not. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)() const>
  void
  depositInterp(EBAMRCellData&              a_meshData,
                const ParticleContainer<P>& a_particles,
                const DepositionType        a_depositionType,
                const bool                  a_forceIrregNGP = false);

  /*!
    @brief Deposit particles on the mesh, keeping the original particle width everywhere and depositing directly
    into every level.
    @param[out] a_meshData       Mesh data. 
    @param[in]  a_particles      Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType Specification of deposition kernel (e.g., CIC)
    @param[in]  a_forceIrregNGP  Force NGP deposition in irregular cells or not. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)() const>
  void
  depositHalo(EBAMRCellData&              a_meshData,
              const ParticleContainer<P>& a_particles,
              const DepositionType        a_depositionType,
              const bool                  a_forceIrregNGP = false);

  /*!
    @brief Deposit particles on the mesh, using an NGP scheme for coarse-grid particles on the refinement boundary. 
    @param[out] a_meshData       Mesh data. 
    @param[in]  a_particles      Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType Specification of deposition kernel (e.g., CIC)
    @param[in]  a_forceIrregNGP  Force NGP deposition in irregular cells or not. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)() const>
  void
  depositHaloNGP(EBAMRCellData&              a_meshData,
                 const ParticleContainer<P>& a_particles,
                 const DepositionType        a_depositionType,
                 const bool                  a_forceIrregNGP = false);

  /*!
    @brief Deposit particles on the mesh.
    @details In this version the particles abutting the refinement boundary from the coarse side are
    deposited with the fine-grid resolution. 
    @param[out] a_meshData       Mesh data. 
    @param[in]  a_particles      Particle container. Must be in "usable state" for deposition.
    @param[in]  a_depositionType Specification of deposition kernel (e.g., CIC)
    @param[in]  a_forceIrregNGP  Force NGP deposition in irregular cells or not. 
  */
  template <class P, class Ret, Ret (P::*MemberFunc)() const>
  void
  depositTransition(EBAMRCellData&              a_meshData,
                    const ParticleContainer<P>& a_particles,
                    const DepositionType        a_depositionType,
                    const bool                  a_forceIrregNGP = false);

  /*!
    @brief Support function for transferring particles that lie in the transition zone around a refinement boundary
    @param[inout] a_particles Input particles
    @param[in] a_depositionType Deposition type.
  */
  template <class P>
  void
  transferMaskParticlesTransition(ParticleContainer<P>& a_particles, const DepositionType a_depositionType) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBAMRParticleMeshImplem.H>

#endif
