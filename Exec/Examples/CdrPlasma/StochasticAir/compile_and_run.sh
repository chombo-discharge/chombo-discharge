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
INPUT="positive2d.inputs Driver.max_steps=20"

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
    mpiexec --report-bindings --bind-to none --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex $INPUT
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'CdrPlasmaTagger::tagCells(EBAMRTags)' \
		       'CdrPlasmaFieldTagger::computeTracers()' \
		       'CdrPlasmaStepper::advanceReactionNetwork(Vector<LD<EBCellFAB>x5, LD<EBCellFAB>, Real, Real, int)' \
		       'CdrPlasmaStepper::computeCdrDiffusionCell(Vector<LD<EBCellFAB>* >x2, LD<EBCellFAB, int, Real)' \
		       'CdrPlasmaStepper::computeCdrDiffusionEb(Vector<LD<BaseIVFAB<Real> >* >x2, LD<BaseIVFAB<Real> >, Real, int)' \
		       'CdrPlasmaStepper::computeCdrDriftVelocities(Vector<LD<EBCellFAB>*>x2, LD<EBCellFAB>, int, Real)' \
		       'CdrPlasmaStepper::computeCdrFluxes(Vector<LD<BaseIVFAB<Real> >*>x6, LD<BaseIVFAB<Real> >, Real)' \
		       'CdrPlasmaStepper::computeCdrDomainFluxes(Vector<LD<DomainFluxIFFAB>*>x6, DomainFluxIFFAB, Real, int)' \
		       'CdrPlasmaStepper::computeElectricField(EBAMRFluxData, phase, EBAMRCellData)' \
		       'CdrPlasmaStepper::extrapolateToDomainFaces(LD<DomainFluxIFFAB>, phase, LD<EBCellFAB>, int)' \
		       'CdrPlasmaStepper::initialSigma()' \
		       'CdrPlasmaStepper::projectFlux(LD<BaseIVFAB<Real> >x2, int)' \
		       'CdrPlasmaStepper::projectDomain(EBAMRIFDatax2)' \
		       'CdrPlasmaStepper::resetDielectricCells(EBAMRIVData)' \
		       'CdrPlasmaStepper::computeElectrodeCurrent()' \
		       'CdrPlasmaStepper::computeDielectricCurrent()' \
		       'CdrPlasmaStepper::computeDomainCurrent()' \
		       'CdrPlasmaStepper::computePhysicsPlotVars' \
		       'CdrPlasmaGodunovStepper::advance(Real)' \
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
