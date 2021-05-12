# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBFastFluxRegister.H\"/\#include \"CD_EbFastFluxRegister.H\"/g' $i
    sed -i 's/\#include \"EBFastFluxRegisterImplem.H\"/\#include \"CD_EbFastFluxRegisterImplem.H\"/g' $i
    
    sed -i 's/\#include <EBFastFluxRegister.H>/\#include <CD_EbFastFluxRegister.H>/g' $i
    sed -i 's/\#include <EBFastFluxRegisterImplem.H>/\#include <CD_EbFastFluxRegisterImplem.H>/g' $i
done

# Find and replace all instances of EBFastFluxRegister. This is now EbCentroidInterpolation
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/EBFastFluxRegister/EbFastFluxRegister/g' $i
done

# Move files
mv src/AmrMesh/EBFastFluxRegister.H    src/AmrMesh/CD_EbFastFluxRegister.H
mv src/AmrMesh/EBFastFluxRegister.cpp  src/AmrMesh/CD_EbFastFluxRegister.cpp
    
