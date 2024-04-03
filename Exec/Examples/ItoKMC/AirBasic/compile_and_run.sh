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
INPUT="example.inputs Driver.max_steps=5 Driver.verbosity=10 ItoKMCGodunovStepper.verbosity=10"
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
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.ex $INPUT
    cp time.table time.table.serial

    # Run OpenMP version
    ./program${DIM}d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex $INPUT
    cp time.table time.table.omp

    # Run MPI version
    mpiexec --report-bindings -n $NCORES --bind-to core ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex $INPUT
    cp time.table.0 time.table.mpi

    # Run MPI+OpenMP version
    mpiexec --report-bindings --bind-to none --map-by slot:PE=$OMP_NUM_THREADS -n $NPROCS ./program${DIM}d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex $INPUT
    cp time.table.0 time.table.hybrid
fi

# Kernel comparison for Source/AmrMesh
if $PROFILE
then
    for PATTERN in 'ParticleOps::getPhysicalParticlesPerCell' \
		       'ParticleOps::getComputationalParticlesPerCell' \
		       'ParticleOps::copy(ParticleContainer<P> x2)' \
		       'ParticleOps::copyDestructive(ParticleContainer<P> x2)' \
		       'ParticleOps::sum(ParticleContainer<P>)' \
		       'ParticleOps::sum(ParticleContainer<P>)' \
		       'ParticleOps::removeParticles(ParticleContainer<P>, std::function<bool(const P&)>)' \
		       'ParticleOps::transferParticles(ParticleContainer<P>, ParticleContainer<P>& std::function<bool(const P&)>)' \
		       'ParticleOps::setData(ParticleContainer<P>, std::function<void(P&)>)' \
		       'ParticleOps::setValue(ParticleContainer<P>, Real' \
		       'ParticleOps::setValue(ParticleContainer<P>, RealVect' \
		       'ParticleContainer<P>::sortParticles()' \
		       'ParticleContainer<P>::organizeParticlesByPatch' \
		       'ParticleContainer<P>::organizeParticlesByPatch' \
		       'ParticleContainer<P>::addParticles(ParticleContainer)' \
		       'ParticleContainer<P>::addParticlesDestructive(ParticleContainer)' \
		       'ParticleContainer<P>::transferParticles(AMRParticles)' \
		       'ParticleContainer::evictInvalidParticles' \
		       'ParticleContainer<P>::remap' \
		       'ParticleContainer<P>::preRegrid' \
		       'ParticleContainer<P>::setValue' \
		       'ParticleContainer<P>::getNumberOfValidParticlesLocal' \
		       'ParticleContainer<P>::getNumberOfOutcastParticlesLocal' \
		       'ParticleContainer<P>::getNumberOfMaskParticlesLocal' \
		       'ParticleContainer<P>::copyMaskParticles(level)' \
		       'ParticleContainer<P>::transferMaskParticles(level)' \
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
