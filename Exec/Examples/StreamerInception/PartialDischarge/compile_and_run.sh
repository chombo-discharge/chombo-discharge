export CH_TIMER=1
export DIM=2
export NCORES=8
export NPROCS=1
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_SCHEDULE="dynamic,4"

COMPILE=true
RUN=true
PROFILE=true

# Compile for serial, OpenMP, flat MPI, and MPI+OpenMP
if $COMPILE
then
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=FALSE DIM=$DIM
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE MPI=FALSE DIM=$DIM
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=TRUE DIM=$DIM
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE MPI=TRUE DIM=$DIM
fi

if $RUN
then
    # Run serial version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.ex example.inputs
    cp time.table time.table.serial
    cp report.txt serial.txt

    # Run OpenMP version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex example.inputs
    cp time.table time.table.omp
    cp report.txt omp.txt

    # # Run MPI version
    mpiexec --report-bindings -n $NCORES --bind-to core ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex example.inputs
    cp time.table.0 time.table.mpi
    cp report.txt mpi.txt    

    # # Run MPI+OpenMP version
    mpiexec --report-bindings --bind-to none --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex example.inputs
    cp time.table.0 time.table.hybrid
    cp report.txt hybrid.txt        
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'StreamerInceptionStepper::seedUniformParticles' \
		       'StreamerInceptionStepper::seedIonizationParticles' \
		       'StreamerInceptionStepper::inceptionIntegrateEuler' \
		       'StreamerInceptionStepper::inceptionIntegrateTrapezoidal' \
		       'StreamerInceptionStepper::townsendTrackEuler' \
		       'StreamerInceptionStepper::townsendTrackTrapezoidal' \
		       'StreamerInceptionStepper::computeRdot' \
		       'StreamerInceptionStepper::rewindTracerParticles' \
		       'StreamerInceptionStepper::resetTracerParticles' \
		       'StreamerInceptionStepper::computeBackgroundIonizationStationary' \
		       'StreamerInceptionStepper::computeDetachmentStationary' \
		       'StreamerInceptionStepper::computeFieldEmissionStationary' \
		       'StreamerInceptionStepper::computeFieldEmission' \
		       'StreamerInceptionStepper::evaluateFunction(level)' \
		       'StreamerInceptionStepper::computeInceptionVoltageVolume' \
		       'StreamerInceptionStepper::computeMinimumInceptionVoltage' \
		       'StreamerInceptionStepper::computeCriticalVolumeStationary' \
		       'StreamerInceptionStepper::computeCriticalVolumeTransient' \
		       'StreamerInceptionStepper::computeCriticalAreaStationary' \
		       'StreamerInceptionStepper::computeCriticalAreaTransient' \
		       'StreamerInceptionStepper::computeIonizationVolumeStationary' \
		       'StreamerInceptionStepper::computeIonizationVolumeTransient' \
		       'StreamerInceptionStepper::computeIonVelocity' \
		       'StreamerInceptionStepper::computeIonDiffusion' \
		       'StreamerInceptionTagger::computeTracerField()' \
		       'StreamerInceptionTagger::tagCells()' \
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
