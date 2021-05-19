for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/electrode/Electrode/g' $i
    sed -i 's/_Electrode/_electrode/g' $i
    sed -i 's/Electrode_/electrode_/g' $i
    sed -i 's/ Electrode/ electrode/g' $i
    sed -i 's/electrode&/Electrode&/g' $i
    sed -i 's/is_live/isLive/g' $i
    sed -i 's/get_function/getImplicitFunction/g' $i
    sed -i 's/get_fraction/getFraction/g' $i
    sed -i 's/electrode(/Electrode(/g' $i
    sed -i 's/class electrode/class Electrode/g' $i
    sed -i 's/electrode::/Electrode::/g' $i
    sed -i 's/Electrodeelectrode/Electrode/g' $i
    sed -i 's/ElectrodeBcs/electrodeBcs/g' $i
    sed -i 's/\.Electrode/\.electrode/g' $i
    sed -i 's/\"Electrode\"/\"electrode\"/g' $i
    sed -i 's/define_Electrode/define_electrode/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"Electrode.H\"/\#include <CD_Electrode.H>/g' $i
    sed -i 's/\#include <Electrode.H>/\#include <CD_Electrode.H>/g' $i

    sed -i 's/\#include \"electrode.H\"/\#include <CD_Electrode.H>/g' $i
    sed -i 's/\#include <electrode.H>/\#include <CD_Electrode.H>/g' $i
done

# Move files
mv Source/Geometry/electrode.H   Source/Geometry/CD_Electrode.H
mv Source/Geometry/electrode.cpp Source/Geometry/CD_Electrode.cpp
