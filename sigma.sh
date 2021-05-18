# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"sigma_solver.H\"/\#include \"CD_SigmaSolver.H\"/g' $i
    sed -i 's/\#include <sigma_solver.H>/\#include <CD_SigmaSolver.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/sigma_solver/SigmaSolver/g' $i
    sed -i 's/get_flux/getFlux/g' $i
    sed -i 's/set_sigma/setSigma/g' $i
    sed -i 's/reset_cells/resetCells/g' $i
    sed -i 's/compute_rhs/computeRHS/g' $i
done

# Move files
mv Source/sigma_solver/sigma_solver.H   Source/sigma_solver/CD_SigmaSolver.H
mv Source/sigma_solver/sigma_solver.cpp Source/sigma_solver/CD_SigmaSolver.cpp
mv Source/sigma_solver Source/SigmaSolver
