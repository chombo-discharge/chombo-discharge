for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/dielectric/Dielectric/g' $i
    sed -i 's/\_Dielectric/\_dielectric/g' $i
    sed -i 's/Dielectric_/dielectric_/g' $i
    sed -i 's/ Dielectric/ dielectric/g' $i
    sed -i 's/dielectric&/Dielectric&/g' $i
    sed -i 's/is_live/isLive/g' $i
    sed -i 's/get_permittivity/getPermittivity/g' $i
    sed -i 's/dielectric(/Dielectric(/g' $i
    sed -i 's/class dielectric/class Dielectric/g' $i
    sed -i 's/dielectric::/Dielectric::/g' $i
    sed -i 's/Dielectricdielectric/Dielectric/g' $i
    sed -i 's/DielectricBcs/dielectricBcs/g' $i
    sed -i 's/\.Dielectric/\.dielectric/g' $i
    sed -i 's/\"Dielectric\"/\"dielectric\"/g' $i
    sed -i 's/define_Dielectric/define_dielectric/g' $i
    sed -i 's/rod_Dielectric/rod_dielectric/g' $i
    sed -i 's/Dielectrics\./dielectrics\./g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"Dielectric.H\"/\#include <CD_Dielectric.H>/g' $i
    sed -i 's/\#include <Dielectric.H>/\#include <CD_Dielectric.H>/g' $i

    sed -i 's/\#include \"dielectric.H\"/\#include <CD_Dielectric.H>/g' $i
    sed -i 's/\#include <dielectric.H>/\#include <CD_Dielectric.H>/g' $i
done

# Move files
mv Source/Geometry/dielectric.H   Source/Geometry/CD_Dielectric.H
mv Source/Geometry/dielectric.cpp Source/Geometry/CD_Dielectric.cpp
