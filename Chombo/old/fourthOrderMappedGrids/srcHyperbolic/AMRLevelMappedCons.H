#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMAPPEDCONS_H_
#define _AMRLEVELMAPPEDCONS_H_

#include "AMRLevelCons.H"
#include "LevelMappedConsOperator.H"
#include "MultiBlockCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"

#include "NamespaceHeader.H"

/// AMR Level for fourth-order solution to conservation law
/**
 */
class AMRLevelMappedCons : public AMRLevelCons
{
  friend class AMRLevelMappedConsFactory ;
public:
  /// Constructor
  /**
   */
  AMRLevelMappedCons();

  /// Destructor
  /**
   */
  virtual ~AMRLevelMappedCons();

  /// Advance by one timestep -- RK4
  /**
   */
  virtual Real advance();

  /// Set a_U at time a_time.
  /**
   */
  virtual void setData(LevelData<FArrayBox>& a_U,
                       Real a_time) const;

  virtual void setDataMapped(LevelData<FArrayBox>& a_U,
                             Real a_time,
                             bool a_includeJ = true) const;

#ifdef CH_USE_HDF5

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

  virtual void getPlotData(LevelData<FArrayBox>& a_plot_data) const;

  virtual int numPlotComps() const;

  /// write mapped-grid info
  virtual void writeMappedPlotFile() const;

#endif

  /// Compute dt using initial data
  virtual Real computeInitialDt();

  void coordinateSystem(MultiBlockCoordSysFactory* a_coordSysFact);

  /// Set whether to use a source term
  void useSourceTerm(bool a_useSourceTerm);

  /// Set the source term
  void sourceTerm(const LevelSourceTerm* const a_sourceTerm);

protected:

  virtual void setDefaultValues();

  // Transfer settings from this to m_levelConsOperator:  called by levelSetup()
  virtual void transferSettingsToLevelOp();

  // Transfer settings from a_amrConsPtr to this:  called by define()
  virtual void transferSettingsFromAMRLevel(AMRLevelMappedCons* a_amrConsPtr);

  /// Compute dt using m_Unew
  virtual Real computeNewDt();

  // Level integrator
  //  LevelMappedConsOperator* m_levelConsOperatorPtr;

  // Pointer to the class defining the physics of the problem
  //  MOLPhysics* m_molPhysics;

  // whether to compute timestep from cell data (the old way)
  bool m_dtFromCells;

  int m_levelStep;

  bool m_useSourceTerm;

  LevelSourceTerm* m_sourceTermPtr;

  // Coordinate system
  //  NewFourthOrderCoordSys* m_coordSysPtr;
  MultiBlockCoordSys* m_coordSysPtr;
  MultiBlockCoordSysFactory* m_coordSysFactPtr;

  MultiBlockLevelGeom m_geom;
  MultiBlockLevelExchangeAverage* m_mblexPtr;

  /// need this in AMRLevel in order to generate mapped-grid output
  std::string m_plotfile_prefix;

  /// whether to write mapped-grid geometry file
  static bool s_writeMap;

  /// whether to write J
  static bool s_writeJ;

  /// whether to write volume
  static bool s_writeVolume;

  /// whether to write error
  static bool s_writeError;
};

#include "NamespaceFooter.H"

#endif
