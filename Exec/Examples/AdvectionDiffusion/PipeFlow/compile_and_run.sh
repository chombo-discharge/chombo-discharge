export CH_TIMER=1
export DIM=2
export NCORES=8
export NPROCS=1
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_SCHEDULE="dynamic,4"

COMPILE=false
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

    # Run OpenMP version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex example.inputs
    cp time.table time.table.omp

    # # Run MPI version
    mpiexec --report-bindings -n $NCORES --bind-to core ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex example.inputs
    cp time.table.0 time.table.mpi

    # # Run MPI+OpenMP version
    mpiexec --report-bindings --bind-to none --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex example.inputs
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'CdrSolver::computeAdvectionFlux(LD<EBFluxFAB>' \
		       'CdrSolver::computeDiffusionFlux(LD<EBFluxFAB>' \
		       'CdrSolver::computeAdvectionDiffusionFlux(EBAMRFluxData' \
		       'CdrSolver::resetDomainFlux(EBAMRFluxData)' \
		       'CdrSolver::fillDomainFlux(LD<EBFluxFAB>, int)' \
		       'CdrSolver::computeDivergenceIrregular(LD<EBCellFAB>' \
		       'CdrSolver::defineInterpolationStencils()' \
		       'CdrSolver::initialDataParticles' \
		       'CdrSolver::hybridDivergence(LD<EBCellFAB>, LD<BaseIVFAB<Real> >, LD<BaseIVFAB<Real> >, int)' \
		       'CdrSolver::interpolateFluxToFaceCentroids(LD<EBFluxFAB>, int)' \
		       'CdrSolver::computeAdvectionDt()' \
		       'CdrSolver::computeDiffusionDt()' \
		       'CdrSolver::computeAdvectionDiffusionDt()' \
		       'CdrSolver::computeSourceDt(Real, Real)' \
		       'CdrSolver::weightedUpwind()' \
		       'CdrSolver::computeMass(EBAMRCellData)' \
		       'CdrSolver::extrapolateAdvectiveFluxToEB(EBAMRIVData)' \
		       'CdrSolver::smoothHeavisideFaces(EBAMRFluxData, EBAMRCellData)' \
		       'CdrSolver::fillGwn(EBAMRFluxData, Real)' \
		       'CdrCTU::computeAdvectionDt()' \
		       'CdrCTU::advectToFaces(EBAMRFluxData, EBAMRCellData, Real)' \
		       'CdrGodunov::computeAdvectionDt()' \
		       'CdrGodunov::advectToFaces(EBAMRFluxDat, EBAMRCellData, Real)' \
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
