/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   Driver.H
  @brief  Declaration of the driver class for chombo-discharge
  @author Robert Marskar
*/

#ifndef CD_Driver_H
#define CD_Driver_H

// Std includes
#include <iostream>
#include <fstream>

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_TimeStepper.H>
#include <CD_AmrMesh.H>
#include <CD_CellTagger.H>
#include <mfis.H>
#include <geo_coarsener.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Main class for time/space advancement of streamer equations. This class runs a timestepper and controls input/output, 
  regrid intervals etc. 
*/
class Driver {
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
	 const RefCountedPtr<AmrMesh>&                a_amr,
	 const RefCountedPtr<CellTagger>&            a_cellTagger = RefCountedPtr<CellTagger>(NULL),
	 const RefCountedPtr<geo_coarsener>&          a_geoCoarsen = RefCountedPtr<geo_coarsener>(NULL));

  /*!
    @brief Constructor
  */
  ~Driver();

  /*!
    @brief Setup and run simulations
    @param[in] a_inputFile Input file containing parameters for the simulation. 
  */
  void setupAndRun(const std::string a_inputFile);
  
  /*!
    @brief Run a simulation for a specific time, or until a maximum number of steps have been performed.
    @param[in] a_startTime Starting time of simulation (seconds)
    @param[in] a_endTime   Simulation end time (seconds)
    @param[in] a_maxSteps  Maximum number of steps.
  */
  void run(const Real a_startTime, const Real a_endTime, const int a_maxSteps);

  /*!
    @brief Setup function
    @param[in] a_inputFile      Input file containing all tunable simulation parameters.
    @param[in] a_initialRegrids Number of initial regrids to perform
    @param[in] a_restart        Do a restart or not
    @param[in] a_restartFile    Name of restart file (if restarting from a previous time step)
  */
  void setup(const std::string a_inputFile, const int a_initialRegrids, const bool a_restart, const std::string a_restartFile);

protected:

  /*!
    @brief Enum for how to print memory to std::cout (or pout)
  */
  enum class MemoryReportMode {
    overall,
    unfreed,
    allocated
  };

  /*!
    @brief Index space
  */
  RefCountedPtr<mfis> m_multifluidIndexSpace;

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
  RefCountedPtr<geo_coarsener> m_geoCoarsen;

  /*!
    @brief Time code
  */
  TimeCode m_timeCode;

  /*!
    @brief Memory report mode
  */
  MemoryReportMode m_memoryReportMode;

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
    @brief Geometric tag depth for dielectrics
  */
  int m_dielectricTagsDepth;

  /*!
    @brief Geometric tag depth for conductors
  */
  int m_conductorTagsDepth;

  /*!
    @brief Geometric tag depth for gas-solid interfaces
  */
  int m_gasSolidInterfaceTagDepth;

  /*!
    @brief Geometry tag depth for solid-solid interfaces
  */
  int m_solidSolidInterfaceTagDepth;

  /*!
    @brief Geometric tag depth for gas-Dielectric interfaces
  */
  int m_gasDielectricInterfaceTagDepth;

  /*!
    @brief Geometric tag depth for gas-conductor interfaces
  */
  int m_gasConductorInterfaceTagDepth;

  /*!
    @brief Maximum number of steps
  */
  int m_maxSteps;

  /*!
    @brief Maximum plot depth
  */
  int m_maxPlotDepth;

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
    @brief Cached tags. Needed when we regrid and change box ownership.
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
    @brief Plot level sets
  */
  bool m_plotLevelset;

  /*!
    @brief Set amr
  */
  void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Sanity check
  */
  void sanityCheck();

  /*!
    @brief Parse options
  */
  void parseOptions();

  /*!
    @brief Parse runtime options
  */
  void parseRuntimeOptions();

  /*!
    @brief Parse geometrical refinement
  */
  void parseGeometryRefinement();

  /*!
    @brief Parse plot variables
  */
  void parsePlotVariables();

  /*!
    @brief Parse option for geometry generation
  */
  void parseGeometryGeneration();

  /*!
    @brief Create output directories
  */
  void createOutputDirectories();

  /*!
    @brief Allocate internal storage for Driver
  */
  void allocateInternals();
  
  /*!
    @brief Cache tags
  */
  void cacheTags(const EBAMRTags& a_tags);

  /*!
    @brief Deallocate internal storage for Driver
  */
  void deallocateInternals();

  /*!
    @brief Write the geometry to file
  */
  void writeEBIS();

  /*!
    @brief Get geometric tags
  */
  void getGeometryTags();

  /*!
    @brief Get boxes and number of cells for everything
  */
  void getLoadsAndBoxes(long long& a_myPoints,
			long long& a_myPointsGhosts,
			long long& a_myBoxes,
			long long& a_totalPoints,
			long long& a_totalPointsGhosts,
			long long& a_totalBoxes,
			Vector<long long>& a_localLevelBoxes,
			Vector<long long>& a_allLevelBoxes,
			Vector<long long>& a_localLevelPoints,
			Vector<long long>& a_totalLevelPoints,
			const int&                       a_finestLevel,
			const Vector<DisjointBoxLayout>& a_grids);

  /*!
    @brief Write a grid report
  */
  void gridReport();

  /*!
    @brief Number formatting game. This just prints the number 123456 as 123,456
  */
  const std::string numberFmt(const long long a_number, char a_sep = ',') const;

  /*!
    @brief Number formatting game. This just prints the number 123456 as 123,456
  */
  const Vector<std::string> numberFmt(const Vector<long long> a_number, char a_sep = ',') const;

  /*!
    @brief Memory report
  */
  void memoryReport(const MemoryReportMode a_mode);

  /*!
    @brief Read checkpoint file
  */
  void readCheckpointFile(const std::string& a_restartFile);

  /*!
    @brief Read vector data up to and including a_step
  */
  void readVectorData(HDF5HeaderData& a_header, Vector<Real>& a_data, const std::string a_name, const int a_elements);

  /*!
    @brief Do a regrid step. This regrid method regrids data inside each solver. The exception is the Poisson solver
    whose solution is recomputed. The same goes for stationary RTE methods. 
  */
  void regrid(const int a_lmin, const int a_lmax, const bool a_useInitialData);

  /*!
    @brief Regrid internal storage for this class
  */
  void regridInternals(const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Write a regrid report
  */
  void regridReport(const Real a_totalTime,
		     const Real a_tagTime,
		     const Real a_baseRegridTime,
		     const Real a_solverRegridTime);
  
  /*!
    @brief Set the computational geometry
  */
  void setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set the time stepper
  */
  void setTimeStepper(const RefCountedPtr<TimeStepper>& a_timeStepper);

  /*!
    @brief Set the cell tagger
  */
  void setCellTagger(const RefCountedPtr<CellTagger>& a_cellTagger);

  /*!
    @brief Set the geo coarsener
  */
  void setGeoCoarsen(const RefCountedPtr<geo_coarsener>& a_geoCoarsen);

  /*!
    @brief Setup function
  */
  void setupFresh(const int a_initialRegrids);

  /*!
    @brief Set for restart
  */
  void setupForRestart(const int a_initialRegrids, const std::string a_restartFile);

  /*!
    @brief Set up for geometry only
  */
  void setupGeometryOnly();

  /*!
    @brief Check if restart file exists
  */
  void checkRestartFile(const std::string a_restartFile) const;

  /*!
    @brief Write a time step report
  */
  void stepReport(const Real a_startTime, const Real a_endTime, const int a_maxSteps);

  /*!
    @brief Write a plot file
  */
  void writePlotFile();

  /*!
    @brief Write a regrid file
  */
  void writeRegridFile();

  /*!
    @brief Write a restart file
  */
  void writeRestartFile();

  /*!
    @brief Write a crash file
  */
  void writeCrashFile();

  /*!
    @brief Write a plot file
  */
  void writePlotFile(const std::string a_filename);

  /*!
    @brief Write the memory usage to file
  */
  void writeMemoryUsage();

  /*!
    @brief Write computational loads
  */
  void writeComputationalLoads();

  /*!
    @brief Write a checkpoint file
  */
  void writeCheckpointFile();
  
  /*!
    @brief Write checkpoint data
  */
  void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write tags
  */
  void writeCheckpointTags(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write checkpoint Realm loads
  */
  void writeCheckpointRealmLoads(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write checkpoint data
  */
  void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Read Realm loads. 
  */
  void readCheckpointRealmLoads(Vector<long int>& a_loads, HDF5Handle& a_handle, const std::string a_realm, const int a_level);

  /*!
    @brief Write vector data to a HDF5Header. This writes up to a_step
  */
  void writeVectorData(HDF5HeaderData& a_header, const Vector<Real>& a_data, const std::string a_name, const int a_elem);

  /*!
    @brief Write the geometry to file. 
  */
  void writeGeometry();

  /*!
    @brief Tag cells for refinement. This computes cell tags and global tags (union of cell tags with geometric tags);
  */
  bool tagCells(Vector<IntVectSet>& a_allTags, EBAMRTags& a_cellTags);

  /*!
    @brief Get number of plot variables
  */
  int getNumberOfPlotVariables() const;

  /*!
    @brief Get names for the plot variables
  */
  Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Write internal plot data
  */
  void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write cell tags to file
  */
  void writeTags(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write MPI ranks to output
  */
  void writeRanks(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write level sets
  */
  void writeLevelset(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Return the finest tag level
  */
  int getFinestTagLevel(const EBAMRTags& a_cellTags) const;

  /*!
    @brief Rebuild Parmparse
  */
  void rebuildParmParse() const;
};

#include <CD_NamespaceFooter.H>

#endif
