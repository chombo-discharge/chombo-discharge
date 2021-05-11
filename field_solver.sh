# Inclusion guards
for i in `find . -name "*.H" -type f`; do
    sed -i 's/\#include \"field_solver.H\"/\#include \"CD_FieldSolver.H\"/g' $i
    sed -i 's/\#include \"field_solver_multigrid.H\"/\#include \"CD_FieldSolverMultigrid.H\"/g' $i
    sed -i 's/\#include \"field_solver_factory.H\"/\#include \"CD_FieldSolverFactory.H\"/g' $i
    sed -i 's/\#include \"field_solver_factoryI.H\"/\#include \"CD_FieldSolverFactoryImplem.H\"/g' $i
done
for i in `find . -name "*.cpp" -type f`; do
    sed -i 's/\#include \"field_solver.H\"/\#include \"CD_FieldSolver.H\"/g' $i
    sed -i 's/\#include \"field_solver_multigrid.H\"/\#include \"CD_FieldSolverMultigrid.H\"/g' $i
    sed -i 's/\#include \"field_solver_factory.H\"/\#include \"CD_FieldSolverFactory.H\"/g' $i
    sed -i 's/\#include \"field_solver_factoryI.H\"/\#include \"CD_FieldSolverFactoryImplem.H\"/g' $i
done

# Rename field_solver_multigrid to FieldSolverMultigrid and field_solver to FieldSolver
for i in `find . -name "*.H" -type f`; do
    sed -i 's/field_solver_multigrid/FieldSolverMultigrid/g' $i
    sed -i 's/field_solver/FieldSolver/g' $i
    sed -i 's/field_solver_factory/FieldSolverFactory/g' $i
    sed -i 's/field_stepper/FieldStepper/g' $i
done
for i in `find . -name "*.cpp" -type f`; do
    sed -i 's/field_solver_multigrid/FieldSolverMultigrid/g' $i
    sed -i 's/field_solver/FieldSolver/g' $i
    sed -i 's/field_stepper/FieldStepper/g' $i
done
for i in `find . -name "*.inputs" -type f`; do
    sed -i 's/field_solver_multigrid/FieldSolverMultigrid/g' $i
done
for i in `find . -name "*.options" -type f`; do
    sed -i 's/field_solver_multigrid/FieldSolverMultigrid/g' $i
done

# Update gnumakefiles everywhere
for i in `find . -name "*GNUmakefile" -type f`; do
    sed -i 's/field_solver/FieldSolver/g' $i
done

# Rename files
mv src/field_solver/field_solver.H src/field_solver/CD_FieldSolver.H
mv src/field_solver/field_solver.cpp src/field_solver/CD_FieldSolver.cpp

mv src/field_solver/field_solver_multigrid.H src/field_solver/CD_FieldSolverMultigrid.H
mv src/field_solver/field_solver_multigrid.cpp src/field_solver/CD_FieldSolverMultigrid.cpp
mv src/field_solver/field_solver_multigrid.options src/field_solver/CD_FieldSolverMultigrid.options

mv src/field_solver/field_solver_factory.H src/field_solver/CD_FieldSolverFactory.H
mv src/field_solver/field_solver_factoryI.H src/field_solver/CD_FieldSolverFactoryImplem.H

# Move field_solver folder to FieldSolver folder
mv src/field_solver src/FieldSolver

