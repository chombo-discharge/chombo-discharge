# This is an example file for testing OpenMP functionality. It should be removed before a pull request.
#
# Specific examples that target folders whose files have been threaded are:
#
# Source/AmrMesh -> Exec/Examples/Electrostatics/MechShaft/compile_and_run.sh
# Source/ConvectionDiffusionReaction -> Exec/Examples/AdvectionDiffusion/PipeFlow
# Source/Driver -> Exec/Examples/CdrPlasma/DeterministicAir/compile_and_run.sh
# Source/Electrostatics -> Exec/Examples/Electrostatics/ProfiledSurface
# Source/Elliptic -> Exec/Examples/Electrostatics/ProfiledSurface
# Source/Geometry -> Exec/Tests/Geometry/RodPlaneProfile
# Source/ItoDiffusion -> Exec/Tests/BrownianWalker/DriftDiffusion
# Source/Particle -> Exec/Examples/ItoKMC/AirBasic
#
# Physics/AdvectionDiffusion -> Exec/Examples/AdvectionDiffusion/PipeFlow
# Physics/BrownianWalker -> Exec/Tests/BrownianWalker/DriftDiffusion
# Physics/ItoKMC -> Exec/Examples/ItoKMC/AirDBD => How to advance KMC in a thread-safe way????
# Physics/StreamerInception -> Exec/Examples/StreamerInception/PartialDis
# Physics/TracerParticle -> Exec/Tests/TracerParticle/CoaxialCable

#
# NOTE: ALL FOLDERS SHOULD BE THREADED BEFORE MERGING.
# NOTE: CHECK ALL CRITICAL PERFORMANCE BOTTLENECKS BEFORE MERGING

export CH_TIMER=1
export DIM=2
export NCORES=8
export NPROCS=1
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_SCHEDULE="dynamic"

COMPILE=true
RUN=true
PROFILE=true
INPUT="regression.inputs"

# Compile for serial, OpenMP, flat MPI, and MPI+OpenMP
if $COMPILE
then
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=FALSE DIM=$DIM
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE  MPI=FALSE DIM=$DIM
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=TRUE  DIM=$DIM
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE  MPI=TRUE  DIM=$DIM
fi

if $RUN
then
    # Run serial version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.ex $INPUT
    cp time.table time.table.serial

    # Run OpenMP version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex $INPUT
    cp time.table time.table.omp

    # # Run MPI version
    mpiexec --report-bindings -n $NCORES --bind-to core ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex $INPUT
    cp time.table.0 time.table.mpi

    # # Run MPI+OpenMP version
    mpiexec --report-bindings --bind-to none --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex $INPUT
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'CdrSolver::computeAdvectionFlux(LD<EBFluxFAB>' \
		   ; do

	if grep -q "${PATTERN}" time.table.serial
	then
	    echo $PATTERN
	    grep -n "${PATTERN}" time.table.serial | head -1
	    grep -n "${PATTERN}" time.table.omp | head -1
	    grep -n "${PATTERN}" time.table.mpi | head -1
	    grep -n "${PATTERN}" time.table.hybrid | head -1
	    echo ""
	fi
    done
fi
