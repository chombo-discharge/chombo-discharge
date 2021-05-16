# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"time_stepper.H\"/\#include \"CD_TimeStepper.H\"/g' $i
    sed -i 's/\#include <time_stepper.H>/\#include <CD_TimeStepper.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/time_stepper/TimeStepper/g' $i
    sed -i 's/_timestepper/_timeStepper/g' $i
    sed -i 's/initial_data/initialData/g' $i
    sed -i 's/time_code/TimeCode/g' $i
    sed -i 's/post_initialize/postInitialize/g' $i
    sed -i 's/post_checkpoint_setup/postCheckpointSetup/g' $i
    sed -i 's/get_redistribution_regsize/getRedistributionRegSize/g' $i
    sed -i 's/write_checkpoint_data/writeCheckpointData/g' $i
    sed -i 's/read_checkpoint_data/readCheckpointData/g' $i
    sed -i 's/get_checkpoint_loads/getCheckpointLoads/g' $i
    sed -i 's/compute_dt/computeDt/g' $i
    sed -i 's/pre_regrid/preRegrid/g' $i
    sed -i 's/need_to_regrid/needToRegrid/g' $i
    sed -i 's/need_to_regrid/needToRegrid/g' $i
    sed -i 's/LoadBalancing_Realm/loadBalanceThisRealm/g' $i
    sed -i 's/LoadBalancing_boxes/loadBalanceBoxes/g' $i
    sed -i 's/post_regrid/postRegrid/g' $i
    sed -i 's/_timecode/_timeCode/g' $i
    sed -i 's/_finest_level/_finestLevel/g' $i
    sed -i 's/_plotvar_names/_plotVariableNames/g' $i
    sed -i 's/TimeSteppers/timesteppers/g' $i
    sed -i 's/TimeCode::advection_diffusion/TimeCode::AdvectionDiffusion/g' $i
    sed -i 's/TimeCode::advection/TimeCode::Advection/g' $i
    sed -i 's/TimeCode::diffusion/TimeCode::Diffusion/g' $i
    sed -i 's/TimeCode::source/TimeCode::Source/g' $i
    sed -i 's/TimeCode::relaxation_time/TimeCode::RelaxationTime/g' $i
    sed -i 's/TimeCode::restricted/TimeCode::Restricted/g' $i
    sed -i 's/TimeCode::hardcap/TimeCode::Hardcap/g' $i
    sed -i 's/TimeCode::error/TimeCode::Error/g' $i
    sed -i 's/TimeCode::physics/TimeCode::Physics/g' $i  
done

# Move files
mv src/Driver/time_stepper.cpp src/Driver/CD_TimeStepper.cpp
mv src/Driver/time_stepper.H   src/Driver/CD_TimeStepper.H

# Update Py and makefiles
for i in `find . -type f \( -iname \*GNUmakefile -o -iname \*.py \)`; do
    sed -i 's/time_stepper/TimeStepper/g' $i
    sed -i 's/TimeSteppers/time_steppers/g' $i # So that we don't fuck with py setup (yet)
done

# enum class TimeCode {
#   Advection,
#   Diffusion,
#   AdvectionDiffusion,  
#   Source,
#   RelaxationTime,
#   Restricted,
#   Hardcap,
#   Error,
#   Physics
# };
