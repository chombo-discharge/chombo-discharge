/*!
  @file   mfis.H
  @brief  Multi-fluid index space
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _mfindexspace_
#define _mfindexspace_

#include "type_definitions.H"

#include <GeometryService.H>
#include <EBIndexSpace.H>
#include <RealVect.H>
#include <MFIndexSpace.H>

/*!
  @brief Two-fluid index space
*/
class mfis {
public:


  /*!
    @brief Weak constructor. Must call define later
  */
  mfis();

  /*!
    @brief Destructor
  */
  ~mfis();

  /*!
    @brief Define function
  */
  virtual void define(const Box                     & a_domain,
		      const RealVect                & a_origin,
		      const Real                    & a_dx,
		      const Vector<GeometryService*>& a_geoservers,
		      int                             a_nCellMax                              = -1,
		      int                             a_maxCoarsenings                        = -1,
		      bool                            a_fixOnlyFirstPhaseRegNextToMultiValued = false);

  /*!
    @brief Get a particular ebis
  */
  const RefCountedPtr<EBIndexSpace>& get_ebis(const phase::which_phase a_phase) const;

  /*!
    @brief Get a particular ebis
  */
  const RefCountedPtr<EBIndexSpace>& get_ebis(const int a_phase) const;

  /*!
    @brief Get a particular ebis
  */
  RefCountedPtr<EBIndexSpace>& get_ebis(const phase::which_phase a_phase);

  /*!
    @brief Get ebis
  */
  RefCountedPtr<EBIndexSpace>& get_ebis(const int a_phase);

  /*!
    @brief Get the interface region. This is defined as the intersection between the irregular cells of the 
    gas phase and the other phases. 
  */
  IntVectSet interface_region(const ProblemDomain& a_domain) const;

  /*!
    @brief Get number of phases
  */
  int num_phases() const;

protected:

  /*!
    @brief EBIndexSpaces
  */
  Vector<RefCountedPtr<EBIndexSpace> > m_ebis;

  /*!
    @brief MFIndexSpace
  */
  RefCountedPtr<MFIndexSpace> m_mfis;
};

#endif
