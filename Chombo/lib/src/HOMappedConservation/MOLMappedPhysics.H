#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MOLMAPPEDPHYSICS_H_
#define _MOLMAPPEDPHYSICS_H_

#include <string>
using std::string;

#include "MOLPhysics.H"
#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

class HDF5HeaderData;

///
/**
   Extends the base class MOLPhysics to include awareness of mapped grids
*/
class MOLMappedPhysics : public MOLPhysics
{
public:
  /// Constructor
  /**
   */
  explicit MOLMappedPhysics();

  /// Destructor
  /**
   */
  virtual ~MOLMappedPhysics();

  /// Copy the data from us into a ptr of the same (base) type
  virtual void copyTo(MOLMappedPhysics* a_newMOLMP) const;

  /// Set the current coordinate system (default implementation - do nothing)
  /**
   */
  virtual void setCurrentCoordSys(const NewCoordSys* a_coordSys);

  /// Set the current time (default implementation - do nothing)
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);


  /// Set face-averaged primitive state on boundary faces
  /**
     Default implementation calls m_bc->primBC.
     Need to call setUnitNormalBasisPtr() before this.
  */
  virtual void primBC(FArrayBox&             a_WGdnv,
                      const FArrayBox&       a_WLeft,
                      const FArrayBox&       a_WRight,
                      const FArrayBox&       a_W,
                      const Real&            a_time,
                      const int&             a_dir);


  virtual void setUnitNormalBasisPtr(const FluxBox *const a_unitNormalBasisPtr);

  virtual void mappedArtVisc(FluxBox&         a_NtF,
                             const FArrayBox& a_U,
                             const FluxBox&   a_N,
                             const FArrayBox& a_J,
                             const Real&      a_alpha,
                             const Real&      a_currentTime,
                             const Box&       a_box);

  virtual void getAux(FArrayBox& a_auxFab,
                      const FArrayBox& a_UFab);

  virtual Vector<string> auxNames();

protected:

  Real m_time;
  NewFourthOrderCoordSys* m_coordSys;

  bool m_haveUnitNormalBasis;
  FluxBox* m_unitNormalBasisPtr;

private:

  // Disallowed for all the usual reasons
  MOLMappedPhysics& operator=(const MOLMappedPhysics& a_input);
  MOLMappedPhysics(const MOLMappedPhysics& a_input);
};

#include "NamespaceFooter.H"
#endif
