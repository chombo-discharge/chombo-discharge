/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   Driver.H
  @brief  Declaration of the simulation driver class for chombo-discharge. 
  @author Robert Marskar
*/

#ifndef CD_Driver_H
#define CD_Driver_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_TimeStepper.H>
#include <CD_AmrMesh.H>
#include <CD_CellTagger.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_GeoCoarsener.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Main class for time/space advancement of streamer equations. 
  @details This class runs a simulation and controls input/output, regrid intervals etc. 
*/
class Driver
{
public:
  /*!
    @brief Weak constructor. Not allowed.
  */
  Driver() = delete;

  /*!
    @brief Full constructor. 
    @param[in] a_computationalGeometry Computational geometry. 
    @param[in] a_TimeStepper          The time stepper
    @param[in] a_amr                   AmrMesh core class. 
    @param[in] a_cellTagger            Cell tagger which flags cell for refinement. 
    @param[in] a_geoCoarsen            geo_coarsen class which can remove tags along the geometry. 
  */
  Driver(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry,
         const RefCountedPtr<TimeStepper>&           a_timeStepper,
         const RefCountedPtr<AmrMesh>&               a_amr,
         const RefCountedPtr<CellTagger>&            a_cellTagger = RefCountedPtr<CellTagger>(nullptr),
         const RefCountedPtr<GeoCoarsener>&          a_geoCoarsen = RefCountedPtr<GeoCoarsener>(nullptr));

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~Driver();

  /*!
    @brief Setup and run simulations
    @param[in] a_inputFile Input file containing parameters for the simulation. 
  */
  void
  setupAndRun(const std::string a_inputFile);

  /*!
    @brief Run a simulation for a specific time, or until a maximum number of steps have been performed.
    @param[in] a_startTime Starting time of simulation
    @param[in] a_endTime   Simulation end time
    @param[in] a_maxSteps  Maximum number of steps.
  */
  void
  run(const Real a_startTime, const Real a_endTime, const int a_maxSteps);

  /*!
    @brief Setup function
    @param[in] a_inputFile      Input file containing all tunable simulation parameters.
    @param[in] a_initialRegrids Number of initial regrids to perform
    @param[in] a_restart        Do a restart or not
    @param[in] a_restartFile    Name of restart file (if restarting from a previous time step)
  */
  void
  setup(const std::string a_inputFile,
        const int         a_initialRegrids,
        const bool        a_restart,
        const std::string a_restartFile);

protected:
  /*!
    @brief Index space
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief AMR
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Geometry
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief Time stepper
  */
  RefCountedPtr<TimeStepper> m_timeStepper;

  /*!
    @brief Cell tagger
  */
  RefCountedPtr<CellTagger> m_cellTagger;

  /*!
    @brief Geometric coarsener
  */
  RefCountedPtr<GeoCoarsener> m_geoCoarsen;

  /*!
    @brief Input file
  */
  std::string m_inputFile;

  /*!
    @brief Name of realm where Driver allocates his data. This is always Realm::Primal.
  */
  std::string m_realm;

  /*!
    @brief Geometry generation method -- can choose between native Chombo or chombo-discharge
  */
  std::string m_geometryGeneration;

  /*!
    @brief Output directory for output files. 
  */
  std::string m_outputDirectory;

  /*!
    @brief Prefix for all output files
  */
  std::string m_outputFileNames;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Irregular tag growth
  */
  int m_irregTagGrowth;

  /*!
    @brief Geometry scan level. 
    @details This is a special option for when we use chombo-discharge's geometry generation method. 
  */
  int m_geoScanLevel;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Regrid interval
  */
  int m_regridInterval;

  /*!
    @brief Checkpoint interval
  */
  int m_checkpointInterval;

  /*!
    @brief Plot interval
  */
  int m_plotInterval;

  /*!
    @brief Geometric tag depth
  */
  int m_geometricTagsDepth;

  /*!
    @brief Geometric tag depth for conductors
  */
  int m_conductorTagsDepth;

  /*!
    @brief Geometric tag depth for dielectrics
  */
  int m_dielectricTagsDepth;

  /*!
    @brief Maximum number of steps
  */
  int m_maxSteps;

  /*!
    @brief Maximum plot depth
  */
  int m_maxPlotLevel;

  /*!
    @brief Maximum plot depth
  */
  int m_maxCheckpointDepth;

  /*!
    @brief Number of initial regrids
  */
  int m_initialRegrids;

  /*!
    @brief Number of ghost cells to plot
  */
  int m_numPlotGhost;

  /*!
    @brief Restart step
  */
  int m_restartStep;

  /*!
    @brief Grow tagged cells
  */
  int m_growTags;

  /*!
    @brief Tags
  */
  Vector<IntVectSet> m_geomTags;

  /*!
    @brief Tags
  */
  EBAMRTags m_tags;

  /*!
    @brief Cached tags. This is done so that when we regrid and change box ownership we can copy the tags around. 
  */
  EBAMRBool m_cachedTags;

  /*!
    @brief Simulation start time
  */
  Real m_startTime;

  /*!
    @brief Simulation stop time
  */
  Real m_stopTime;

  /*!
    @brief Time increment
  */
  Real m_dt;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief 
  */
  Real m_capacitance;

  /*!
    @brief MPI timer
  */
  Real m_wallClockStart;

  /*!
    @brief MPI timer
  */
  Real m_wallClockOne;

  /*!
    @brief MPI timer
  */
  Real m_wallClockTwo;

  /*!
    @brief Angle refinement threshold
  */
  Real m_refineAngle;

  /*!
    @brief Profile or not
  */
  bool m_profile;

  /*!
    @brief Turn on/off geometry coarsening below the base level
  */
  bool m_doCoarsening;

  /*!
    @brief Special option for when geometric tags are changed during a simulation. 
    @details This becomes = true in parseGeometryRefinement every time the geometric refinement criteria changed, and it always becomes false after a regrid. 
  */
  bool m_needsNewGeometricTags;

  /*!
    @brief Write regrid files or not
  */
  bool m_writeRegridFiles;

  /*!
    @brief Write restart files or not
  */
  bool m_writeRestartFiles;

  /*!
    @brief Allow grid coarsening
  */
  bool m_allowCoarsening;

  /*!
    @brief Do init load balancing.
  */
  bool m_doInitLoadBalancing;

  /*!
    @brief Write memory report
  */
  bool m_writeMemory;

  /*!
    @brief Write loads to file or not
  */
  bool m_writeLoads;

  /*!
    @brief Restart or not
  */
  bool m_restart;

  /*!
    @brief Geometry only
  */
  bool m_geometryOnly;

  /*!
    @brief Use memory as load balance for EBIS
  */
  bool m_ebisMemoryLoadBalance;

  /*!
    @brief Plot tagged cells
  */
  bool m_plotTags;

  /*!
    @brief Plot MPI ranks
  */
  bool m_plotRanks;

  /*!
    @brief Plot patch loads
  */
  bool m_plotLoads;

  /*!
    @brief Plot level sets
  */
  bool m_plotLevelset;

  /*!
    @brief Set amr
    @param[in] a_amrMesh AmrMesh
  */
  void
  setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh);

  /*!
    @brief Sanity check
  */
  void
  sanityCheck();

  /*!
    @brief Parse options
    @details Called at the beginning of the simulation. 
  */
  void
  parseOptions();

  /*!
    @brief Parse runtime options
    @details Called every time step. 
  */
  void
  parseRuntimeOptions();

  /*!
    @brief Parse geometrical refinement 
    @details This parses settings for how deep to refine geometries, either based on curvature or a hard limit for electrodes/dielectrics. 
  */
  void
  parseGeometryRefinement();

  /*!
    @brief Parse plot variables
  */
  void
  parsePlotVariables();

  /*!
    @brief Parse option for geometry generation. 
    @details This sets the geometry-generation load balancing method to either use Chombo or chombo-discharge. 
  */
  void
  parseGeometryGeneration();

  /*!
    @brief Parse irregular tag growth
    @details This parses m_irregTagGrowth, which is used to grow cut-cell tags. 
  */
  void
  parseIrregTagGrowth();

  /*!
    @brief Create output directories
  */
  void
  createOutputDirectories();

  /*!
    @brief Allocate internal storage for Driver. 
  */
  void
  allocateInternals();

  /*!
    @brief Cache tags. This is used
  */
  void
  cacheTags(const EBAMRTags& a_tags);

  /*!
    @brief Write the geometry to file
  */
  void
  writeEBIS();

  /*!
    @brief Get geometric tags
    @details This fills m_geomTags with irregular cell tags, using information that was passed into Driver from the input script. This includes
    adding tags that refine surfaces to a certain depth, or tags where the curvature of the surface exceeds a specified value. 
  */
  void
  getGeometryTags();

  /*!
    @brief Get statics about a DisjointBoxLaout. This will compute the number of grid cells and patches both locally and globally. 
    @param[out] a_numLocalCells       Local number of grid cells, excluding ghost cells. 
    @param[out] a_numLocalCellsGhosts Local number of grid cells, including ghost cells. 
    @param[out] a_numLocalBox         Local number of grid patches (over all levels)
    @param[out] a_numTotalCells       Total number of grid cells, excluding ghost cells. 
    @param[out] a_numTotalCellsGhosts Total number of grid cells, including ghost cells. 
    @param[out] a_numTotalBox         Total number of grid patches
    @param[out] a_numLocalLevelBoxes  Local number of grid boxes (per level)
    @param[out] a_numTotalLevelBoxes  Total number of grid boxes (per level)
    @param[out] a_numLocalLevelCells  Local number of grid cells (per level)
    @param[out] a_numTotalLevelCells  Total number of grid cells (per level)
    @param[in]  a_finestLevel         Finest AMR level
    @param[in]  a_grids               Input grids
  */
  void
  getCellsAndBoxes(long long&                       a_numLocalCells,
                   long long&                       a_numLocalCellsGhosts,
                   long long&                       a_numLocalBoxes,
                   long long&                       a_numTotalCells,
                   long long&                       a_numTotalCellsGhosts,
                   long long&                       a_numTotalBoxes,
                   Vector<long long>&               a_numLocalLevelBoxes,
                   Vector<long long>&               a_numTotalLevelBoxes,
                   Vector<long long>&               a_numLocalLevelCells,
                   Vector<long long>&               a_numTotalLevelCells,
                   const int&                       a_finestLevel,
                   const Vector<DisjointBoxLayout>& a_grids);

  /*!
    @brief Print regrid report. 
  */
  void
  gridReport();

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpoint file
    @details This will read grids and computational loads from an HDF5 file and restart the simulation from that. 
    @param[in] a_restartFile Which file to read
  */
  void
  readCheckpointFile(const std::string& a_restartFile);
#endif

  /*!
    @brief Do a regrid step. 
    @param[in] a_lmin           Coarsest level that is allowed to change during regrid
    @param[in] a_lmax           Finest level that can be generated during the regrid step. 
    @param[in] a_useInitialData If true, fill solvers with initial data after making the new grids. 
  */
  void
  regrid(const int a_lmin, const int a_lmax, const bool a_useInitialData);

  /*!
    @brief Regrid internal storage for this class
    @param[in] a_oldFinestLevel Finest level before the regrid
    @param[in] a_newFinestLevel Finest level before the regrid
  */
  void
  regridInternals(const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Set the computational geometry
    @param[in] a_computationalGeometry Computational geometry
  */
  void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set the time stepper
    @param[in] a_timeStepper Time stepper
  */
  void
  setTimeStepper(const RefCountedPtr<TimeStepper>& a_timeStepper);

  /*!
    @brief Set the cell tagger
    @param[in] a_cellTagger Cell tagger
  */
  void
  setCellTagger(const RefCountedPtr<CellTagger>& a_cellTagger);

  /*!
    @brief Set the geo coarsener
    @param[in] a_geoCoarsen Coarsener for cut-cell tags. 
  */
  void
  setGeoCoarsen(const RefCountedPtr<GeoCoarsener>& a_geoCoarsen);

  /*!
    @brief Setup function
  */
  void
  setupFresh(const int a_initialRegrids);

  /*!
    @brief Set for restart
  */
#ifdef CH_USE_HDF5
  void
  setupForRestart(const int a_initialRegrids, const std::string a_restartFile);
#endif

  /*!
    @brief Set up for geometry only
  */
  void
  setupGeometryOnly();

  /*!
    @brief Check if restart file exists
  */
  void
  checkRestartFile(const std::string a_restartFile) const;

  /*!
    @brief Write a time step report
  */
  void
  stepReport(const Real a_startTime, const Real a_endTime, const int a_maxSteps);

  /*!
    @brief Write a plot file
    @param[in] a_file File name
    @details This can write an arbitrary folder. E.g. filename = /crash/abc.hdf5
  */
  void
  writePlotFile(const std::string a_filename);

  /*!
    @brief Write a plot file. This writes to plt/
  */
  void
  writePlotFile();

  /*!
    @brief Write a regrid file.
    @details This writes a regular plot file to /regrid
  */
  void
  writePreRegridFile();

  /*!
    @brief Write a regrid file.
    @details This writes a regular plot file to /regrid
  */
  void
  writePostRegridFile();

  /*!
    @brief Write a restart file
    @details This writes a regular plot file to /restart
  */
  void
  writeRestartFile();

  /*!
    @brief Write a crash file
    @details This writes a regular plot file to /crash
  */
  void
  writeCrashFile();

  /*!
    @brief Write the memory usage to file
  */
  void
  writeMemoryUsage();

  /*!
    @brief Write computational loads
  */
  void
  writeComputationalLoads();

  /*!
    @brief Write a checkpoint file
  */
  void
  writeCheckpointFile();

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data. 
    @param[in] a_handle HDF5 file
    @param[in] a_level  Grid level
    @details This will call writeCheckpointTags and writeCheckpointRealmLoads
  */
  void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write tags to checkpoint file
    @param[in] a_handle HDF5 file
    @param[in] a_level  Grid level
  */
  void
  writeCheckpointTags(HDF5Handle& a_handle, const int a_level);
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write computational loads to checkpoint file
    @param[in] a_handle HDF5 file
    @param[in] a_level  Grid level
  */
  void
  writeCheckpointRealmLoads(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Write checkpoint data
  */
#ifdef CH_USE_HDF5
  void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Read Realm loads. 
  */
#ifdef CH_USE_HDF5
  void
  readCheckpointRealmLoads(Vector<long int>& a_loads,
                           HDF5Handle&       a_handle,
                           const std::string a_realm,
                           const int         a_level);
#endif

  /*!
    @brief Write the geometry to file. 
  */
  void
  writeGeometry();

  /*!
    @brief Tag cells for refinement. This computes cell tags and global tags (union of cell tags with geometric tags);
  */
  bool
  tagCells(Vector<IntVectSet>& a_allTags, EBAMRTags& a_cellTags);

  /*!
    @brief Get number of plot variables
    @return Returns the number of internal plot variables that Driver will write to file. 
  */
  int
  getNumberOfPlotVariables() const;

  /*!
    @brief Get the string identifiers for the internal plot variables from Driver
  */
  Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Write internal plot data
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Starting component
    @param[in]    a_level  AMR level
    @details This will copy Driver's plot data into a_output, starting on component a_comp. 
  */
  virtual void
  writePlotData(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept;

  /*!
    @brief Write cell tags to file
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Starting component
    @param[in]    a_level  AMR level
    @details This will write cell tags data into a_output, starting on component a_comp. 
  */
  virtual void
  writeTags(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept;

  /*!
    @brief Write MPI ranks to output
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Starting component
    @param[in]    a_level  AMR level
    @details This adds data to a_output where the data in each grid patch is the MPI rank. 
  */
  virtual void
  writeRanks(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept;

  /*!
    @brief Write level sets
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Starting component
    @details This adds data to a_output where the cell-centered data is the level-set value. 
    @note This writes directly to a_output.
  */
  virtual void
  writeLevelset(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept;

  /*!
    @brief Write computational loads to output
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Starting component
    @param[in]    a_level  AMR level
    @details This adds data to a_output where the data in each grid patch is the MPI rank. 
  */
  virtual void
  writeLoads(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept;

  /*!
    @brief Return the finest tag level
  */
  int
  getFinestTagLevel(const EBAMRTags& a_cellTags) const;

  /*!
    @brief Rebuild Parmparse
  */
  void
  rebuildParmParse() const;
};

#include <CD_NamespaceFooter.H>

#endif
