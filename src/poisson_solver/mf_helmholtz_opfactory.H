/*!
  @file mf_helmholtz_opfactory.H
  @brief Factory class for creating multifluid helmholtz operators
  @author Robert Marskar
*/

#ifndef _MF_HELMHOLTZ_OPFACTORY_
#define _MF_HELMHOLTZ_OPFACTORY_

#include "type_definitions.H"
#include "mf_helmholtz_op.H"

#include <AMRMultiGrid.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <ProblemDomain.H>
#include <RefCountedPtr.H>
#include <EBIndexSpace.H>
#include <MFCellFAB.H>

class mf_helmholtz_opfactory : public AMRLevelOpFactory<LevelData<MFCellFAB> > {
public:

  /*!
    @brief Factory class for creating multifluid helmholtz operators
  */
  mf_helmholtz_opfactory(const RefCountedPtr<EBIndexSpace>&        a_ebis_gas,
			 const RefCountedPtr<EBIndexSpace>&        a_ebis_solid,
			 const Vector<int>&                        a_ref_ratio,
			 const ProblemDomain&                      a_coarsest_domain,
			 const Real&                               a_coarsest_dx,
			 const RealVect&                           a_origin,
			 const RefCountedPtr<BaseDomainBCFactory>& a_domainbc_factory,
			 const RefCountedPtr<BaseEBBCFactory>&     a_ebbc_factory,
			 const IntVect&                            a_ghost_phi,
			 const IntVect&                            a_ghost_rhs,
			 int                                       a_num_levels = -1);

  /*!
    @brief Destructor
  */
  ~mf_helmholtz_opfactory();

  /*!
    @brief Set when to drop to the bottom solver
  */
  static void set_bottom_drop(const int a_bottom_drop);

    /*!
    @brief Set when to drop to the bottom solver
  */
  static void set_relax_type(const int a_relax_type);

  /*!
    @brief Set jump condition
    @param[in] a_gas_perm Gas-side permittivity
    @param[in] a_solid_perm Solid permittivity on jump cells
    @param[in] a_jump Surface charge (divided by eps0)
  */
  virtual void set_jump(const Real a_gas_perm, const EBAMRIVData& a_solid_perm, const EBAMRIVData& a_jump);

  /*!
    @brief Factory method
  */
  virtual MGLevelOp<LevelData<MFCellFAB> >* MGnewOp(const ProblemDomain& a_fine_ebis, int a_depth, bool a_homo_only = true);

  /*!
    @brief Factory method
  */
  virtual AMRLevelOp<LevelData<MFCellFAB> >* AMRnewOp(const ProblemDomain& a_fineindexspace);

  /*!
    @brief Factory method
  */
  virtual mf_helmholtz_op* createOperator(const DisjointBoxLayout&       a_dilboMGLevel,
					  const DisjointBoxLayout&       a_dilboCoarMG,
					  const ProblemDomain&           a_domainMGLevel,
					  const bool&                    a_hasMGObjects,
					  const bool&                    a_layoutChanged,
					  const RealVect&                a_dxMGLevel,
					  const RealVect&                a_dxCoar,
					  const int&                     a_whichLevel,
					  const int&                     a_mgLevel);

  /*!
    @brief Reclaim
  */
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  /*!
    @brief Reclaim
  */
  virtual void AMRreclaim(mf_helmholtz_op* a_reclaim);

  /*!
    @brief Refinement ratio between this level and the coarser level. Returns 1 when there are no coarser levels
  */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

protected:

  /*!
    @brief Gas-side ebis
  */
  RefCountedPtr<EBIndexSpace> m_ebis_gas;

  /*!
    @brief EBIS inside of solids
  */
  RefCountedPtr<EBIndexSpace> m_ebis_solid;
  
  /*!
    @brief Drop to bottom solver
  */
  static int s_drop_bottom;

  /*!
    @brief Relaxation type
  */
  static int s_relax_type;

  /*!
    @brief Defined or not
  */
  bool m_defined;
private:
  
  /*!
    @brief Disallowed weak constructor
  */
  mf_helmholtz_opfactory(){
    MayDay::Error("mf_helmholtz_opfactory:: invalid constructor");
  }
};
#endif
