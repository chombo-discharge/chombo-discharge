# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"mfconductivityop.H\"/\#include <CD_MfConductivityOp.H>/g' $i
    sed -i 's/\#include \"mfconductivityopfactory.H\"/\#include <CD_MfConductivityOpFactory.H>/g' $i
    sed -i 's/\#include <mfconductivityop.H>/\#include <CD_MfConductivityOp.H>/g' $i
    sed -i 's/\#include <mfconductivityopfactory.H>/\#include <CD_MfConductivityOpFactory.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/mfconductivityopfactory/MfConductivityOpFactory/g' $i
    sed -i 's/mfconductivityop/MfConductivityOp/g' $i
done

# Move files
mv Source/Elliptic/mfconductivityopfactory.H   Source/Elliptic/CD_MfConductivityOpFactory.H
mv Source/Elliptic/mfconductivityopfactory.cpp Source/Elliptic/CD_MfConductivityOpFactory.cpp
mv Source/Elliptic/mfconductivityop.H          Source/Elliptic/CD_MfConductivityOp.H
mv Source/Elliptic/mfconductivityop.cpp        Source/Elliptic/CD_MfConductivityOp.cpp
