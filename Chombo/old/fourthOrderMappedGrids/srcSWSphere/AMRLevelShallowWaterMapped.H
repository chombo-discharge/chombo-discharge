#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELSHALLOWWATERMAPPED_H_
#define _AMRLEVELSHALLOWWATERMAPPED_H_

#include "AMRLevelMappedCons.H"

#include "NamespaceHeader.H"

/// AMR Level for fourth-order solution to conservation law
/**
 */
class AMRLevelShallowWaterMapped : public AMRLevelMappedCons
{
  friend class AMRLevelShallowWaterMappedFactory ;
public:
  /// Constructor
  /**
   */
  AMRLevelShallowWaterMapped();

  /// Destructor
  /**
   */
  virtual ~AMRLevelShallowWaterMapped();

  /// Set a_U at time a_time.
  /**
   */
  virtual void setDataMapped(LevelData<FArrayBox>& a_U,
                             Real a_time,
                             bool a_includeJ = true) const;
#ifdef CH_USE_HDF5

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  virtual void getPlotData(LevelData<FArrayBox>& a_plot_data) const;

  virtual int numPlotComps() const;

#endif

protected:

  /*
  /// Compute dt using m_Unew
  virtual Real computeNewDt();
  */

  /// Index within primitive variables for tagging cells
  virtual int indexForTagging();

  virtual void transferSettingsToLevelOp();
};

#include "NamespaceFooter.H"
#endif
