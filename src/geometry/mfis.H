/*!
  @file mfindexspace.H
  @brief Multi-fluid index space
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _mfindexspace_
#define _mfindexspace_

#include "type_definitions.H"

#include <GeometryService.H>
#include <EBIndexSpace.H>
#include <RealVect.H>

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
  void define(const Box                     & a_domain,
	      const RealVect                & a_origin,
	      const Real                    & a_dx,
	      const Vector<GeometryService*>& a_geoservers,
	      int                             a_nCellMax                              = -1,
	      int                             a_maxCoarsenings                        = -1,
	      bool                            a_fixOnlyFirstPhaseRegNextToMultiValued = false);

  /*!
    @brief Get a particular ebis
  */
  const RefCountedPtr<EBIndexSpace>& get_ebis(Phase::WhichPhase a_whichEBIS) const;

  /*!
    @brief Get the interface region. This is defined as the intersection between the irregular cells of the 
    gas phase and the other phases. 
  */
  const IntVectSet interface_region(const ProblemDomain& a_domain) const;

protected:

  /*!
    @brief EBIndexSpaces
  */
  Vector<RefCountedPtr<EBIndexSpace> > m_ebis;
};

#endif
