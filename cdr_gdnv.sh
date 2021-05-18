# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_gdnv.H\"/\#include <CD_CdrGodunov.H>/g' $i
    sed -i 's/\#include <cdr_gdnv.H>/\#include <CD_CdrGodunov.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_solver/CdrSolver/g' $i

done

# Move files
mv src/CdrSolver/cdr_gdnv.H   src/CdrSolver/CD_CdrGodunov.H
mv src/CdrSolver/cdr_gdnv.cpp   src/CdrSolver/CD_CdrGodunov.cpp

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_gdnv/CdrGodunov/g' $i
    sed -i 's/advect_to_faces/advectToFaces/g' $i
    sed -i 's/parse_slopelim/parseSlopeLimiter/g' $i
done
