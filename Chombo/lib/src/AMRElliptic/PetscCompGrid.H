#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCCOMPGRID_H_
#define _PETSCCOMPGRID_H_
#include "Stencil.H"
#include "BCFunc.H"
#include "DisjointBoxLayout.H"
#include "BoxIterator.H"
#include "FourthOrderInterpStencil.H"
#ifdef CH_USE_PETSC
#include <petsc.h>
#endif
#include "NamespaceHeader.H"

typedef enum {GHOST=-1,FINE_COVERED=-2,DIRI_BC=-3,
              NEUM_BC=-4,ANY_BC=-5,UNKNOWN=-6} GID_type; 

std::ostream& operator<< (std::ostream& os, GID_type);

//! \class PetscCompGrid
//! This base class organizes the construction of a PETSc matrix, with an AMR hierarchy
class PetscCompGrid
{

public: 
  //! Destructor.
  virtual ~PetscCompGrid();

#ifdef CH_USE_PETSC

  //! Base class constructor. Called by all subclass constructors.
  PetscCompGrid(int a_dof) : m_gid0(0),m_mat(0),m_Pmat(0),m_CFStencilRad(2),m_writeMatlab(false),
                             m_averageFineSolutionToCoarse(true),m_verbose(0),m_dof(a_dof),
                             m_num_extra_nnz(0),m_repartition(PETSC_FALSE),
                             m_patch_size(0),m_from_petscscat(0)
  {
#if defined(PETSC_USE_LOG)
    PetscLogEventRegister("createMatrix",PETSC_VIEWER_CLASSID,&m_event0);
    PetscLogEventRegister("PetscCompGrid",PETSC_VIEWER_CLASSID,&m_event1);
    PetscLogEventRegister("Prepartsition",PETSC_VIEWER_CLASSID,&m_event2);
#endif
  }
  
  virtual void define( const ProblemDomain &a_cdomain,
                       Vector<DisjointBoxLayout> &a_grids, 
                       Vector<int> &a_refratios, 
                       BCHolder a_bc,
                       const RealVect &a_cdx,
                       int a_numLevels=-1, int a_ibase=0);
  
  virtual void clean();

  Mat getMatrix() const { return m_mat; }
  Mat getPMatrix() const { return m_Pmat; }
  void setMatlab(bool b=true) {m_writeMatlab = b;}
  void setRepartition(bool b=true) {m_repartition = (b ? PETSC_TRUE : PETSC_FALSE);}
  void setVerbose(int a_v){ m_verbose=a_v;}
  void setAverageFineSolutionToCoarse(bool b=true) {m_averageFineSolutionToCoarse = b;}
  
  // main Mat creation routines
  PetscErrorCode createMatrix(int a_makePmat=0);
  // solver uses these, not used in base class, used in app code
  PetscErrorCode putChomboInPetsc( const Vector<LevelData<FArrayBox> * > &rhs, Vec b )const;
  PetscErrorCode putPetscInChombo( Vec b, Vector<LevelData<FArrayBox> * > &rhs )const;
  // derived type needs to tell me how many ghosts it needs
  virtual IntVect getGhostVect()const = 0;
  virtual void addExtraCovered(GID_type,int,const DataIndex&,BaseFab<PetscInt>&) {} // default no op
protected:
  // core stencil ops - base class has basic implementations except for equations
  virtual void createOpStencil(IntVect,int,const DataIndex&,StencilTensor &) = 0;
  virtual void applyBCs(IntVect,int,const DataIndex&,Box,StencilTensor &);
  virtual void InterpToFine(IntVect,int,const DataIndex&,StencilTensor &);
  virtual void InterpToCoarse(IntVect,int,const DataIndex&,StencilTensor &);
  // helper for FourthOrderInterpStencil
  IntVect getCFStencil(const ProblemDomain &a_cdom, const IntVect a_ivc);
  // at end add stencil to matrix - done
  PetscErrorCode AddStencilToMat(IntVect,int,const DataIndex&,StencilTensor &, Mat);
  // utils
  void NodeDefine(StencilNode &a_node, IntVect a_iv, int a_lev, Real a_val)
  {
    a_node.first.setIV(a_iv);
    a_node.first.setLevel(a_lev);
    a_node.second.define(1);
    a_node.second.setValue(a_val); // this sets diagonal
  }
  void setCFCoverMaps(int a_nlev);
  void setCoverMaps(int a_nlev);
  PetscErrorCode permuteDataAndMaps(Vector<StencilTensor> &patchStencil);
  // data
  Vector<ProblemDomain> m_domains;
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_refRatios;
  Vector<RealVect> m_dxs;
public:
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt> > > > m_GIDs;
protected:
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt> > > > m_crsSupportGIDs;
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt> > > > m_fineCoverGIDs;
public:
  PetscInt m_gid0,m_nlocrealpatches,m_patchid0;
protected:
  Mat m_mat;
  Mat m_Pmat;
  BCHolder m_bc;
  // all possible stencils, on box(-m_CFStencilRad:+m_CFStencilRad)^SpaceDim
  int m_CFStencilRad;
  BaseFab<FourthOrderInterpStencil*> m_FCStencils;
  bool m_writeMatlab;
  bool m_averageFineSolutionToCoarse;
  int m_verbose;
  const int m_dof;
#if defined(PETSC_USE_LOG)
  PetscLogEvent m_event0;
  PetscLogEvent m_event1;
  PetscLogEvent m_event2;
#endif
public:
  int m_num_extra_nnz;
  // repartitioning data
  PetscBool m_repartition;
private:
  PetscInt m_patch_size;
  VecScatter m_from_petscscat;
  Vec m_origvec; // buffer to hold Chombo data in a PETSc vector
private:
  // Disallowed for all the usual reasons
  void operator=(const PetscCompGrid& a_input)
  {
    MayDay::Error("invalid operator");
  }
  // Disallowed for all the usual reasons
  PetscCompGrid(const PetscCompGrid& a_input):m_dof(0)
  {
    MayDay::Error("invalid operator");
  }
#endif // ifdef petsc

};

class CompBC : public BCFunction
{
public:
  CompBC() {;}
  CompBC(int a_nSource, IntVect a_nGhosts);
  virtual ~CompBC();
  void define(int a_nSource, IntVect a_nGhosts);

  virtual void createCoefs() = 0;
  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           bool                 a_homogeneous) = 0;

  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           const DataIndex&     a_index,
                           bool                 a_homogeneous) 
  {
    operator()(a_state, a_valid, a_domain, a_dx, a_homogeneous);
  }
 
  IntVect nGhosts()const{return m_nGhosts;}
  int nSources()const{return m_nSources;}
#ifdef CH_USE_PETSC
  PetscReal getCoef(int a_iSrc, int a_iGhost=0);
#else
  Real getCoef(int a_iSrc, int a_iGhost=0);
#endif
protected:
#ifdef CH_USE_PETSC
  PetscReal *m_Rcoefs;
#else
  Vector<Real> m_Rcoefs;
#endif
  IntVect    m_nGhosts;
  int        m_nSources; // degreee is m_nSources-1
  bool       m_isDefined;

};

class ConstDiriBC : public CompBC
{
public:
  ConstDiriBC(int a_nSource=1, IntVect a_nGhosts=IntVect::Unit) : CompBC(a_nSource,a_nGhosts) {}
  virtual void createCoefs();
  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           bool                 a_homogeneous);
};

#include "NamespaceFooter.H"

#endif
