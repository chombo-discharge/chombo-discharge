# #Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"geo_coarsener.H\"/\#include <CD_GeoCoarsener.H>/g' $i
#     sed -i 's/\#include <geo_coarsener.H>/\#include <CD_GeoCoarsener.H>/g' $i
# done

# Move files
# mv Source/Geometry/geo_coarsener.H   Source/Geometry/CD_GeoCoarsener.H
# mv Source/Geometry/geo_coarsener.cpp Source/Geometry/CD_GeoCoarsener.cpp
mv Source/Geometry/geo_coarsener.options Source/Geometry/CD_GeoCoarsener.options

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/geo_coarsener/GeoCoarsener/g' $i
#     sed -i 's/coarsen_tags/coarsenTags/g' $i
#     sed -i 's/get_coarsen_boxes/getCoarsenBoxes/g' $i
#     sed -i 's/get_coarsen_levels/getCoarsenLevels/g' $i
# done
