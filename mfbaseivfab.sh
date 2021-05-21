#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"MFFastFluxReg.H\"/\#include <CD_MFFluxReg.H>/g' $i
    sed -i 's/\#include <MFFastFluxReg.H>/\#include <CD_MFFluxReg.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/get_fastfr_ptr/getFluxRegPointer/g' $i
    sed -i 's/get_fastfr/getFluxReg/g' $i
    sed -i 's/MFFastFluxReg/MFFluxReg/g' $i
done

# Move files
mv Source/Utilities/MFFastFluxReg.H    Source/Utilities/CD_MFFluxReg.H
