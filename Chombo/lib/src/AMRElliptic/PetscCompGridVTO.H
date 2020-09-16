#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCCOMPGRIDVTO_H_
#define _PETSCCOMPGRIDVTO_H_
#include "PetscCompGrid.H"
#include "NamespaceHeader.H"

//! \class PetscCompGridVTO
//! This base class organizes the construction of a PETSc matrix, 
//! and solve with an AMR 
class PetscCompGridVTO : public PetscCompGrid
{
public:

  //! Destructor.
  virtual ~PetscCompGridVTO()
  {   
#ifdef CH_USE_PETSC
    clean();
#endif
  }

#ifdef CH_USE_PETSC
  //! Base class constructor. Called by all subclass constructors.
  PetscCompGridVTO(Real a_al=0., Real a_beta=1.0) : 
    PetscCompGrid(CH_SPACEDIM), m_alpha(a_al), m_beta(a_beta)
  {
  }
  virtual void clean();

  void defineCoefs( Real a_alpha, Real a_beta, 
                    Vector<RefCountedPtr<LevelData<FluxBox> > > &a_mu,
                    Vector<RefCountedPtr<LevelData<FluxBox> > > &a_lambda,
                    Vector<RefCountedPtr<LevelData<FArrayBox> > > &a_a)
  {
    m_alpha = a_alpha;
    m_beta = a_beta;
    m_a = a_a;
    m_eta = a_mu;
    m_lamb = a_lambda;
  }

  virtual IntVect getGhostVect()const
  {
    return IntVect::Unit;
  }
  Real getAlpha()const{return m_alpha;}
  Real getBeta()const{return m_beta;}
  void setAlpha(Real a_b) {m_alpha = a_b;}
  void setBeta(Real a_b) {m_beta = a_b;}
protected:
  virtual void createOpStencil(IntVect,int,const DataIndex&,StencilTensor&);
  virtual void applyBCs(IntVect,int,const DataIndex&,Box,StencilTensor &);
  Real              m_alpha;
  Real              m_beta;
  //coefficient storage. 
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_eta;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_lamb;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_a;
#endif // ifdef petsc
};

class CompGridVTOBC : public CompBC
{
public:
  CompGridVTOBC(int a_nSource=1, IntVect a_nGhosts=IntVect::Unit) : CompBC(a_nSource,a_nGhosts)
  {
    for (int comp=0; comp<SpaceDim; comp++)
      {
        m_bcDiri[0][0][comp] = true;
        m_bcDiri[0][1][comp] = true;
        m_bcDiri[1][0][comp] = true;
        m_bcDiri[1][1][comp] = true;
      }
  }
  // destructor
  virtual ~CompGridVTOBC() 
  {
  }


  void setDiri(int dir, int side, int comp, bool a_diri = true) {m_bcDiri[dir][side][comp] = a_diri;}
  bool isDiri(int dir, int side, int comp) { return m_bcDiri[dir][side][comp]; }
  virtual void createCoefs();
  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           bool                 a_homogeneous);

  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           const DataIndex&     a_index,
                           bool                 a_homogeneous)
  {
    operator()(a_state, a_valid, a_domain, a_dx, a_homogeneous);
  }


  bool m_bcDiri[CH_SPACEDIM][2][CH_SPACEDIM]; // dim, side, Comp: false = Neumann, true: Diri
};

#include "NamespaceFooter.H"
#endif
