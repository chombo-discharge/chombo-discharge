export NCORES=12
export CH_TIMER=1
export OMP_NUM_THREADS=$NCORES
export OMP_PLACES=cores
export OMP_SCHEDULE="static,1"

# Compile for serial, OpenMP, flat MPI, and MPI+OpenMP
# make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=FALSE
# make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE MPI=FALSE
# make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=FALSE MPI=TRUE
# make -s -j$NCORES OPT=HIGH DEBUG=FALSE OPENMPCC=TRUE MPI=TRUE


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
mpiexec -n 1 --bind-to socket ./program3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex example.inputs FieldSolverMultigrid.gmg_exit_tol=1.2
cp time.table.0 time.table.hybrid

PATTERN='EBGhostCellInterpolator::define'
grep -n "${PATTERN}" time.table.serial | head -1
grep -n "${PATTERN}" time.table.omp | head -1
grep -n "${PATTERN}" time.table.mpi | head -1
grep -n "${PATTERN}" time.table.hybrid | head -1
echo ""

PATTERN='EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator'
grep -n "${PATTERN}" time.table.serial | head -1
grep -n "${PATTERN}" time.table.omp | head -1
grep -n "${PATTERN}" time.table.mpi | head -1
grep -n "${PATTERN}" time.table.hybrid | head -1
echo ""

PATTERN='EBCoarAve::define'
grep -n "${PATTERN}" time.table.serial | head -1
grep -n "${PATTERN}" time.table.omp | head -1
grep -n "${PATTERN}" time.table.mpi | head -1
grep -n "${PATTERN}" time.table.hybrid | head -1
echo ""

PATTERN='EBFluxRedistribution::define'
grep -n "${PATTERN}" time.table.serial | head -1
grep -n "${PATTERN}" time.table.omp | head -1
grep -n "${PATTERN}" time.table.mpi | head -1
grep -n "${PATTERN}" time.table.hybrid | head -1
echo ""

PATTERN='EBReflux::define'
grep -n "${PATTERN}" time.table.serial | head -1
grep -n "${PATTERN}" time.table.omp | head -1
grep -n "${PATTERN}" time.table.mpi | head -1
grep -n "${PATTERN}" time.table.hybrid | head -1
echo ""
