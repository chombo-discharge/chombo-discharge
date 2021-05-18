# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_bc.H\"/\#include \"CD_CdrBc.H\"/g' $i
    sed -i 's/\#include <cdr_bc.H>/\#include <CD_CdrBc.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_solver/CdrSolver/g' $i

done

# Move files
mv src/CdrSolver/cdr_bc.H   src/CdrSolver/CD_CdrBc.H

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cdr_bc::which_bc/CdrBc/g' $i
    sed -i 's/cdr_bc/CdrBc/g' $i
    sed -i 's/CdrBc::external/CdrBc::External/g' $i
    sed -i 's/CdrBc::wall/CdrBc::Wall/g' $i
    sed -i 's/CdrBc::outflow/CdrBc::Outflow/g' $i
    sed -i 's/CdrBc::extrap/CdrBc::Extrap/g' $i
done
