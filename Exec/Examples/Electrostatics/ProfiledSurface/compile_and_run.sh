export CH_TIMER=1
export DIM=2
export NCORES=12
export NPROCS=1
export OMP_NUM_THREADS=12
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_SCHEDULE="dynamic"

COMPILE=true
RUN=true
PROFILE=true
INPUT="example2d.inputs Driver.max_steps=0"
# Driver.initial_regrids=1 Driver.write_memory=true Driver.write_loads=true FieldStepper.load_balance=true"

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
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.ex $INPUT Driver.output_names=serial
    cp time.table time.table.serial

    # Run OpenMP version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex $INPUT Driver.output_names=openmp
    cp time.table time.table.omp

    # Run MPI version
    mpiexec --report-bindings -n $NCORES --bind-to core ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex $INPUT Driver.output_names=mpi
    cp time.table.0 time.table.mpi

    # Run MPI+OpenMP version
    mpiexec --report-bindings --bind-to none --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex $INPUT Driver.output_names=hybrid
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'FieldSolver::computeDisplacementField(MFAMRCellData, MFAMRCellData)' \
		       'FieldSolver::setPermittivities()' \
		       'FieldSolver::writeMultifluidData' \
		       'FieldSolver::computeLoads(DisjointBoxLayout, int)' \
		       'MFHelmholtzElectrostaticEBBC::defineSinglePhase()' \
		       'EBHelmholtzOp::dotProduct' \
		       'EBHelmholtzOp::defineStencils()' \
		       'EBHelmholtzOp::norm' \
		       'EBHelmholtzOp::refluxFreeAMROperator' \
		       'EBHelmholtzOp::applyOp(LD<EBCellFAB>' \
		       'EBHelmholtzOp::diagonalScale(LD<EBCellFAB>, bool)' \
		       'EBHelmholtzOp::relaxPointJacobi(LD<EBCellFAB>, LD<EBCellFAB>, int)' \
		       'EBHelmholtzOp::relaxGSRedBlack(LD<EBCellFAB>, LD<EBCellFAB>, int)' \
		       'EBHelmholtzOp::relaxGSMultiColor(LD<EBCellFAB>, LD<EBCellFAB>, int)' \
		       'EBHelmholtzOp::computeDiagWeight()' \
		       'EBHelmholtzOp::computeRelaxationCoefficient()' \
		       'EBHelmholtzOp::makeAggStencil()' \
		       'EBHelmholtzOp::computeFlux(LD<EBCellFAB>)' \
		       'DataOps::setValue(LD<MFCellFAB>, Real)' \
		       'DataOps::averageCellToFace(LD<EBFluxFAB, LD<EBCellFAB>, ....' \
		       'DataOps::axby' \
		       'DataOps::incr(LD<EBCellFAB)' \
		       'DataOps::getMaxMin(Real, Real, LD<EBFluxFAB>, int>)' \
		       'DataOps::kappaScale(LD<EBCellFAB>)' \
		       'DataOps::kappaScale(LD<MFCellFAB>)' \
		       'DataOps::scale(LD<MFCellFAB>)' \
		       'DataOps::scale(LD<EBCellFAB>)' \
		       'DataOps::setCoveredValue(LD<EBFluxFAB>, int, Real)' \
		       'DataOps::setValue(LD<MFCellFAB>, Real)' \
		       'DataOps::incr(LD<MFCellFAB)' \
		       'MFHelmholtzOp::applyOp' \
		       'MFHelmholtzOp::relaxPointJacobi' \
		       'MFHelmholtzOp::interpolateCF' \
		       'MFHelmholtzOp::relaxGSRedBlack' \
		       'MFHelmholtzOp::relaxGSMultiColor' \
		       'MFHelmholtzJumpBC::defineStencils()' \
		       'MFHelmholtzJumpBC::buildAverageStencils()' \
		       'MFHelmholtzJumpBC::defineIterators()' \
		       'MFHelmholtzJumpBC::resetBC()' \
		       'MFHelmholtzJumpBC::matchBC(LD<MFCellFAB>, LD<BaseIVFAB<Real>, bool)' \
		       'EBHelmholtzDirichletEBBC::define()' \
		       'MFHelmholtzEBBC::defineMultiPhase()' \
		       'EBHelmholtzOpFactory::defineMultigridLevels()' \
		       'MFHelmholtzDirichletEBBC::defineSinglePhase()' \
		       'MFHelmholtzRobinEBBC::define()' \
		       'EBHelmholtzNeumannEBBC::define()' \
		       'MFHelmholtzOpFactory::defineJump()' \
		       'MFHelmholtzOpFactory::setJump(EBAMRIVData, Real)' \
		       'EBHelmholtzRobinEBBC::define()' \
		       'MFHelmholtzOp::relax' \
		       'AMRMultiGrid::solveNo-InitResid' \
		       'FieldSolverMultigrid::setupSolver()' \
		       'AMRMultiGrid::computeAMRResidual' \
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
