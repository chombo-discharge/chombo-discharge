# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"ebconductivityop.H\"/\#include <CD_EbConductivityOp.H>/g' $i
    sed -i 's/\#include \"ebconductivityopfactory.H\"/\#include <CD_EbConductivityOpFactory.H>/g' $i
    sed -i 's/\#include <ebconductivityop.H>/\#include <CD_EbConductivityOp.H>/g' $i
    sed -i 's/\#include <ebconductivityopfactory.H>/\#include <CD_EbConductivityOpFactory.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/ebconductivityopfactory/EbConductivityOpFactory/g' $i
    sed -i 's/ebconductivityop/EbConductivityOp/g' $i
done

# Move files
mv Source/Elliptic/ebconductivityopfactory.H   Source/Elliptic/CD_EbConductivityOpFactory.H
mv Source/Elliptic/ebconductivityopfactory.cpp Source/Elliptic/CD_EbConductivityOpFactory.cpp
mv Source/Elliptic/ebconductivityop.H          Source/Elliptic/CD_EbConductivityOp.H
mv Source/Elliptic/ebconductivityop.cpp        Source/Elliptic/CD_EbConductivityOp.cpp
