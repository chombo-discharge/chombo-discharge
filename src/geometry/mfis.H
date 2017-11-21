/*!
  @file mfindexspace.H
  @brief Multi-fluid index space
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _mfindexspace_
#define _mfindexspace_

#include <GeometryService.H>
#include <EBIndexSpace.H>
#include <RealVect.H>

/*!
  @brief Global stuff for phases
*/
namespace Phase {
  const int num_phases = 2;
  
  enum WhichPhase {
    Gas = 0,
    Solid = 1
  };
};

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
  const RefCountedPtr<EBIndexSpace> get_ebis(Phase::WhichPhase a_whichEBIS) const;

protected:

  /*!
    @brief EBIndexSpaces
  */
  Vector<RefCountedPtr<EBIndexSpace> > m_ebis;
};

#endif
