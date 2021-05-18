# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_fhdF_F.H\"/\#include <CD_CdrFhdF_F.H>/g' $i
    sed -i 's/\#include \"cdr_musclF_F.H\"/\#include <CD_CdrMusclF_F.H>/g' $i
    sed -i 's/\#include <cdr_fhdF_F.H>/\#include <CD_CdrFhdF_F.H>/g' $i
    sed -i 's/\#include <cdr_musclF_F.H>/\#include <CD_CdrMusclF_F.H>/g' $i
done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/cdr_solver/CdrSolver/g' $i

# done

# # Move files
mv src/CdrSolver/cdr_fhdF_F.H    src/CdrSolver/CD_CdrFhdF_F.H
mv src/CdrSolver/cdr_fhdF.ChF    src/CdrSolver/CD_CdrFhdF.ChF

mv src/CdrSolver/cdr_musclF_F.H    src/CdrSolver/CD_CdrMusclF_F.H
mv src/CdrSolver/cdr_musclF.ChF    src/CdrSolver/CD_CdrMusclF.ChF
