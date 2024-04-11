export CH_TIMER=1
export DIM=2
export NCORES=8
export NPROCS=1
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_SCHEDULE="dynamic"

COMPILE=false
RUN=true
PROFILE=true
INPUT="example.inputs Driver.max_steps=20"

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
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.ex $INPUT >& pout.serial
    cp time.table time.table.serial

    # Run OpenMP version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex $INPUT >& pout.openmp
    cp time.table time.table.omp

    # # Run MPI version
    mpiexec --report-bindings -n $NCORES --bind-to core ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex $INPUT
    cp time.table.0 time.table.mpi

    # # Run MPI+OpenMP version
    mpiexec --report-bindings --bind-to core --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex $INPUT
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'ItoKMCTagger::tagCells'\
		       'ItoKMCStepper::initialSigma' \
		       'ItoKMCStepper::writeNumberOfParticlesPerPatch' \
		       'ItoKMCStepper::computeMobilities(mobilities, E, level, time)' \
		       'ItoKMCStepper::computeDiffusionCoefficients(Vector<LD<EBCellFAB>*>, LD<EBCellFAB>, int, Real)' \
		       'ItoKMCStepper::computeReactiveItoParticlesPerCell(LD<EBCellFAB>, int)' \
		       'ItoKMCStepper::computeReactiveCdrParticlesPerCell(LD<EBCellFAB>, int)' \
		       'ItoKMCStepper::computeReactiveMeanEnergiesPerCell(LD<EBCellFAB>, int)' \
		       'ItoKMCStepper::advanceReactionNetwork(LD<EBCellFAB>x3, int, Real)' \
		       'ItoKMCStepper::reconcileParticles(LevelData<EBCellFAB>x3, int)' \
		       'ItoKMCStepper::reconcilePhotoionization()' \
		       'ItoKMCStepper::reconcileCdrDensities(EBAMRCellDatax2, Real)' \
		       'ItoKMCStepper::reconcileCdrDensities(LD<EBCellFAB>x2, int, Real)' \
		       'ItoKMCStepper::fillSecondaryEmissionEB(full)' \
		       'ItoKMCStepper::computePhysicsDt(LD<EBCellFAB>, LD<EBCellFAB>, int)' \
		       'ItoKMCStepper::loadBalanceParticleRealm(...)' \
		       'ItoKMCStepper::computeEdotJSource(a_dt)' \
		       'ItoKMCStepper::computePhysicsPlotVariables' \
		       'ItoKMCGodunovStepper::advance' \
		       'ItoKMCGodunovStepper::setOldPositions' \
		       'ItoKMCGodunovStepper::copyConductivityParticles' \
		       'ItoKMCGodunovStepper::diffuseParticlesEulerMaruyama' \
		       'ItoKMCGodunovStepper::stepEulerMaruyamaParticles' \
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
