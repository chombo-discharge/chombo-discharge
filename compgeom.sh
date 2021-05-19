#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"computational_geometry.H\"/\#include <CD_ComputationalGeometry.H>/g' $i
    sed -i 's/\#include <computational_geometry.H>/\#include <CD_ComputationalGeometry.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/computational_geometry/ComputationalGeometry/g' $i
    sed -i 's/get_dielectrics/getDielectrics/g' $i
    sed -i 's/get_electrodes/getElectrodes/g' $i
    sed -i 's/set_dielectrics/setDielectrics/g' $i
    sed -i 's/set_electrodes/setElectrodes/g' $i
    sed -i 's/get_eps0/getGasPermittivity/g' $i
    sed -i 's/set_eps0/setGasPermittivity/g' $i
    sed -i 's/get_mfis/getMfIndexSpace/g' $i
    sed -i 's/get_gas_if/getGasImplicitFunction/g' $i
    sed -i 's/get_sol_if/getSolidImplicitFunction/g' $i
    sed -i 's/build_geometries/buildGeometries/g' $i
    sed -i 's/build_geo_from_files/buildGeometriesFromFiles/g' $i
    sed -i 's/build_gas_geoserv/buildGasGeoServer/g' $i
    sed -i 's/build_sol_geoserv/buildSolidGeoServer/g' $i
done

# Move files
mv Source/Geometry/computational_geometry.H   Source/Geometry/CD_ComputationalGeometry.H
mv Source/Geometry/computational_geometry.cpp Source/Geometry/CD_ComputationalGeometry.cpp
