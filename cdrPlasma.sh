 #Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cdr_plasma_tagger.H\"/\#include <CD_CdrPlasmaTagger.H>/g' $i
    sed -i 's/\#include <cdr_plasma_tagger.H</\#include <CD_CdrPlasmaTagger.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/cdr_plasma_tagger/CdrPlasmaTagger/g' $i
    sed -i 's/compute_tracers/computeTracers/g' $i
    sed -i 's/get_tracer_fields/getTracerFields/g' $i
    sed -i 's/refine_cells_box/refineCellsBox/g' $i
    sed -i 's/coarsen_cells_box/coarsenCellsBox/g' $i
    sed -i 's/coarsen_cell/coarsenCell/g' $i
    sed -i 's/refine_cell/refineCell/g' $i
done

# Move files
mv Physics/CdrPlasma/cdr_plasma_tagger.H   Physics/CdrPlasma/CD_CdrPlasmaTagger.H
mv Physics/CdrPlasma/cdr_plasma_tagger.cpp Physics/CdrPlasma/CD_CdrPlasmaTagger.cpp


# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/using namespace physics::cdr_plasma/using namespace Physics::CdrPlasma/g' $i
#     sed -i 's/namespace physics/namespace Physics/g' $i
#     sed -i 's/namespace cdr_plasma/namespace CdrPlasma/g' $i
# done
