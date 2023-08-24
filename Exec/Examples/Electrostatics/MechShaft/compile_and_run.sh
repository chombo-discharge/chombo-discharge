export CH_TIMER=1
export NCORES=12
export NRANKS=1
export OMP_NUM_THREADS=12
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_SCHEDULE="static, 1"

COMPILE=true
RUN=true
PROFILE_AMRMESH_COAR=true
PROFILE_AMRMESH_FINE=false

# Compile for serial, OpenMP, flat MPI, and MPI+OpenMP
if $COMPILE
then
    # make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=FALSE
    # make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE MPI=FALSE
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=TRUE
    make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE MPI=TRUE
fi

if $RUN
then
    # Run serial version
    ./program3d.Linux.64.g++.gfortran.OPTHIGH.ex example.inputs FieldSolverMultigrid.gmg_exit_tol=1.2
    cp time.table time.table.serial

    # Run OpenMP version
    ./program3d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex example.inputs FieldSolverMultigrid.gmg_exit_tol=1.2 FieldSolverMultigrid.gmg_exit_tol=1.2
    cp time.table time.table.omp

    # # Run MPI version
    mpiexec -n $NCORES ./program3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex example.inputs FieldSolverMultigrid.gmg_exit_tol=1.2
    cp time.table.0 time.table.mpi

    # # Run MPI+OpenMP version
    mpiexec -n $NRANKS --bind-to socket ./program3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex example.inputs FieldSolverMultigrid.gmg_exit_tol=1.2
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh folder
# Kernel comparison for Source/AmrMesh folder
if $PROFILE_AMRMESH_COAR
then
    for PATTERN in 'AmrMesh::regridAmr' \
		       'AmrMesh::regridOperators' \
		       'Realm::regridOperators' \
		       'Realm::regridBase' \
		       'PhaseRealm::regridOperators' \
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

if $PROFILE_AMRMESH_FINE
then
    for PATTERN in 'AmrMesh::removeCoveredParticlesIF' \
		       'AmrMesh::removeCoveredParticlesDiscrete' \
		       'AmrMesh::removeCoveredParticlesVoxels' \
		       'AmrMesh::transferCoveredParticlesIF' \
		       'AmrMesh::transferCoveredParticlesDiscrete' \
		       'AmrMesh::transferCoveredParticlesVoxels' \
		       'AmrMesh::transferIrregularParticles' \
		       'AmrMesh::intersectParticlesRaycastIF' \
		       'AmrMesh::intersectParticlesBisectIF' \
		       'AmrMesh::allocate(EBAMRIVData, string, phase::which_phase, int, int)' \
		       'AmrMesh::reallocate(EBAMRIVData, phase::which_phase, int)' \
		       'EBCoarAve::define' \
		       'EBCoarAve::defineCellStencils' \
		       'EBCoarAve::defineFaceStencils' \
		       'EBCoarAve::defineEBStencils' \
		       'EBCoarAve::averageData(LD<EBCellFAB>)' \
		       'EBCoarAve::averageData(LD<EBFluxFAB>)' \
		       'EBCoarAve::averageData(LD<BaseIVFAB>)' \
		       'EBCoarseToFineInterp::define' \
		       'EBCoarseToFineInterp::defineWeights' \
		       'EBCoarseToFineInterp::interpolate(LD<EBCellFAB>)' \
		       'EBCoarseToFineInterp::interpolate(LD<BaseIVFAB<Real>>)' \
		       'EBFluxRedistribution::defineStencils' \
		       'EBFluxRedistribution::defineValidCells' \
		       'EBFluxRedistribution::defineInterfaceCells' \
		       'EBFluxRedistribution::redistributeCoar' \
		       'EBFluxRedistribution::redistributeLevel' \
		       'EBFluxRedistribution::redistributeFine' \
		       'EBGhostCellInterpolator::defineGhostRegions' \
		       'EBGhostCellInterpolator::interpolate(LD<EBCellFAB>' \
		       'EBGradient::computeLevelGradient' \
		       'EBGradient::computeAMRGradient' \
		       'EBGradient::computeNormalDerivative' \
		       'EBGradient::defineLevelStencils' \
		       'EBGradient::defineMasks' \
		       'EBGradient::defineIteratorsEBCF' \
		       'EBGradient::defineStencilsEBCF' \
		       'EBGradient::makeAggStencils' \
		       'EBLeastSquaresMultigridInterpolator::coarseFineInterp' \
		       'EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LD<EBCellFAB>, Interval)' \
		       'EBLeastSquaresMultigridInterpolator::defineGhostRegions' \
		       'EBLeastSquaresMultigridInterpolator::defineCoarseInterp' \
		       'EBLeastSquaresMultigridInterpolator::defineStencilsEBCF' \
		       'EBLeastSquaresMultigridInterpolator::makeAggStencils' \
		       'EBLeastSquaresMultigridInterpolator::regularCoarseFineInterp' \
		       'EBMGProlong::define' \
		       'EBMGProlong::prolongResidual' \
		       'EBMGRestrict::define' \
		       'EBMGRestrict::restrictResidual' \
		       'EBReflux::defineStencils' \
		       'EBReflux::coarsenFluxes' \
		       'EBReflux::refluxIntoCoarse' \
		       'EBReflux::defineRegionsCF' \
		       'IrregAmrStencil::apply(LD<EBCellFABx2, int)' \
		       'IrregAmrStencil::apply(LD<EBCellFAB, int)' \
		       'IrregAmrStencil::apply(LD<BaseIVFAB>, LD<EBCellFAB>, int)' \
		       'IrregStencil::define' \
		       'PhaseRealm::defineVofIterator' \
		       'PhaseRealm::defineLevelSet' \
		       'Realm::defineHaloMasks' \
		       'Realm::defineValidCells' \
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
