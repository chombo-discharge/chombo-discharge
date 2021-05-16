# # Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"driver.H\"/\#include \"CD_Driver.H\"/g' $i
#     sed -i 's/\#include <driver.H>/\#include <CD_Driver.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/driver/Driver/g' $i
#     sed -i 's/set_amr/setAmr/g' $i
#     sed -i 's/setup_and_run/setupAndRun/g' $i
#     sed -i 's/parse_plot_vars/parsePlotVariables/g' $i
#     sed -i 's/parse_geometry_generation/parseGeometryGeneration/g' $i
#     sed -i 's/parse_geo_refinement/parseGeometryRefinement/g' $i
#     sed -i 's/create_output_directories/createOutputDirectories/g' $i
#     sed -i 's/allocate_internals/allocateInternals/g' $i
#     sed -i 's/cache_tags/cacheTags/g' $i
#     sed -i 's/deallocate_internals/deallocateInternals/g' $i
#     sed -i 's/write_ebis/writeEBIS/g' $i
#     sed -i 's/get_geom_tags/getGeometryTags/g' $i
#     sed -i 's/get_loads_and_boxes/getLoadsAndBoxes/g' $i
#     sed -i 's/grid_report/gridReport/g' $i
#     sed -i 's/number_fmt/numberFmt/g' $i
#     sed -i 's/memory_report_mode/MemoryReportMode/g' $i
#     sed -i 's/memory_report/memoryReport/g' $i
#     sed -i 's/read_checkpoint_file/readCheckpointFile/g' $i
#     sed -i 's/read_vector_data/readVectorData/g' $i
#     sed -i 's/regrid_internals/regridInternals/g' $i
#     sed -i 's/regrid_report/regridReport/g' $i
#     sed -i 's/set_computational_geometry/setComputationalGeometry/g' $i
#     sed -i 's/set_time_stepper/setTimeStepper/g' $i
#     sed -i 's/set_cell_tagger/setCellTagger/g' $i
#     sed -i 's/set_geo_coarsen/setGeoCoarsen/g' $i
#     sed -i 's/setup_fresh/setupFresh/g' $i
#     sed -i 's/setup_for_restart/setupForRestart/g' $i
#     sed -i 's/setup_geometry_only/setupGeometryOnly/g' $i
#     sed -i 's/check_restart_file/checkRestartFile/g' $i
#     sed -i 's/step_report/stepReport/g' $i
#     sed -i 's/write_plot_file/writePlotFile/g' $i
#     sed -i 's/write_regrid_file/writeRegridFile/g' $i
#     sed -i 's/write_restart_file/writeRestartFile/g' $i
#     sed -i 's/write_crash_file/writeCrashFile/g' $i
#     sed -i 's/write_memory_usage/writeMemoryUsage/g' $i
#     sed -i 's/write_computational_loads/writeComputationalLoads/g' $i
#     sed -i 's/write_checkpoint_file/writeCheckpointFile/g' $i
#     sed -i 's/write_checkpoint_level/writeCheckpointLevel/g' $i
#     sed -i 's/write_checkpoint_tags/writeCheckpointTags/g' $i
#     sed -i 's/write_checkpoint_Realm_loads/writeCheckpointRealmLoads/g' $i
#     sed -i 's/read_checkpoint_level/readCheckpointLevel/g' $i
#     sed -i 's/read_checkpoint_Realm_loads/readCheckpointRealmLoads/g' $i
#     sed -i 's/write_vector_data/writeVectorData/g' $i
#     sed -i 's/write_geometry/writeGeometry/g' $i
#     sed -i 's/tag_cells/tagCells/g' $i
#     sed -i 's/get_num_plot_vars/getNumberOfPlotVariables/g' $i
#     sed -i 's/get_plot_var_names/getPlotVariableNames/g' $i
#     sed -i 's/write_plot_data/writePlotData/g' $i
#     sed -i 's/write_tags/writeTags/g' $i
#     sed -i 's/write_ranks/writeRanks/g' $i
#     sed -i 's/write_levelset/writeLevelset/g' $i
#     sed -i 's/_get_finest_tag_level/_getFinestTagLevel/g' $i
#     sed -i 's/_input_file/_inputFile/g' $i
#     sed -i 's/_end_time/_endTime/g' $i
#     sed -i 's/_max_steps/_maxSteps/g' $i
#     sed -i 's/_init_regrids/_initialRegrids/g' $i
#     sed -i 's/_restart_file/_restartFile/g' $i
#     sed -i 's/_compgeom/_computationalGeometry/g' $i
#     sed -i 's/_timestepper/_timeStepper/g' $i
#     sed -i 's/_celltagger/_cellTagger/g' $i
#     sed -i 's/_geocoarsen/_geoCoarsen/g' $i
#     sed -i 's/_memory_mode/_memoryReportMode/g' $i
#     sed -i 's/_geometry_generation/_geometryGeneration/g' $i
#     sed -i 's/_output_dir/_outputDirectory/g' $i
#     sed -i 's/_output_names/_outputFileNames/g' $i
#     sed -i 's/_geo_scan_level/_geoScanLevel/g' $i
#     sed -i 's/_regrid_interval/_regridInterval/g' $i
#     sed -i 's/_chk_interval/_checkpointInterval/g' $i
#     sed -i 's/_plot_interval/_plotInterval/g' $i
#     sed -i 's/_geom_tag_depth/_geometricTagsDepth/g' $i
#     sed -i 's/_dielectric_tag_depth/_dielectricTagsDepth/g' $i
#     sed -i 's/_conductor_tag_depth/_conductorTagsDepth/g' $i
#     sed -i 's/_conductor_tag_depth/_conductorTagsDepth/g' $i
#     sed -i 's/_gas_solid_interface_tag_depth/_gasSolidInterfaceTagDepth/g' $i
#     sed -i 's/_solid_solid_interface_tag_depth/_solidSolidInterfaceTagDepth/g' $i
#     sed -i 's/_gas_dielectric_interface_tag_depth/_gasDielectricInterfaceTagDepth/g' $i
#     sed -i 's/_gas_conductor_interface_tag_depth/_gasConductorInterfaceTagDepth/g' $i
#     sed -i 's/_max_steps/_maxSteps/g' $i
#     sed -i 's/_max_plot_depth/_maxPlotDepth/g' $i
#     sed -i 's/_max_chk_depth/_maxCheckpointDepth/g' $i
#     sed -i 's/_init_regrids/_initialRegrids/g' $i
#     sed -i 's/_num_plot_ghost/_numPlotGhost/g' $i
#     sed -i 's/_restart_step/_restartStep/g' $i
#     sed -i 's/_grow_tags/_growTags/g' $i
#     sed -i 's/_geom_tags/_geomTags/g' $i
#     sed -i 's/_cached_tags/_cachedTags/g' $i
#     sed -i 's/_start_time/_startTime/g' $i
#     sed -i 's/_stop_time/_stopTime/g' $i
#     sed -i 's/_wallclock_start/_wallClockStart/g' $i
#     sed -i 's/_wallclock1/_wallClockOne/g' $i
#     sed -i 's/_wallclock2/_wallClockTwo/g' $i
#     sed -i 's/_write_regrid_files/_writeRegridPlotFiles/g' $i
#     sed -i 's/_write_restart_files/_writeRestartPlotFiles/g' $i
#     sed -i 's/_allow_coarsen/_allowCoarsening/g' $i
#     sed -i 's/_write_memory/_writeMemory/g' $i
#     sed -i 's/_write_loads/_writeLoads/g' $i
#     sed -i 's/_geometry_only/_geometryOnly/g' $i
#     sed -i 's/_ebis_memory_LoadBalancing/_ebisMemoryLoadBalance/g' $i
#     sed -i 's/_plot_tags/_plotTags/g' $i
#     sed -i 's/_plot_ranks/_plotRanks/g' $i
#     sed -i 's/_plot_levelset/_plotLevelset/g' $i
#     sed -i 's/_my_level_boxes/_localLevelBoxes/g' $i
#     sed -i 's/_total_level_boxes/_allLevelBoxes/g' $i
#     sed -i 's/_my_level_points/_localLevelPoints/g' $i
#     sed -i 's/_total_level_points/_totalLevelPoints/g' $i
#     sed -i 's/_use_initial_data/_useInitialData/g' $i
#     sed -i 's/_old_finest_level/_oldFinestLevel/g' $i
#     sed -i 's/_new_finest_level/_newFinestLevel/g' $i
#     sed -i 's/_total_time/_totalTime/g' $i
#     sed -i 's/_tag_time/_tagTime/g' $i
#     sed -i 's/_base_regrid_time/_baseRegridTime/g' $i
#     sed -i 's/_solver_regrid_time/_solverRegridTime/g' $i
#     sed -i 's/_all_tags/_allTags/g' $i
#     sed -i 's/_cell_tags/_cellTags/g' $i
# done

# # Move files
# mv src/driver/driver.cpp     src/driver/CD_Driver.cpp
# mv src/driver/driver.H       src/driver/CD_Driver.H
# mv src/driver/driver.options src/driver/CD_Driver.options

# # Move folder
# mv src/driver src/Driver

# # Update Py and makefiles
# for i in `find . -type f \( -iname \*GNUmakefile -o -iname \*.py \)`; do
#     sed -i 's/driver/Driver/g' $i
# done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/get_finest_tag_level/getFinestTagLevel/' $i
done
