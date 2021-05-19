# #Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"geometry_stepper.H\"/\#include <CD_GeometryStepper.H>/g' $i
#     sed -i 's/\#include <geometry_stepper.H>/\#include <CD_GeometryStepper.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/geometry_stepper/GeometryStepper/g' $i
# done

# # Move files
# mv Physics/geometry/geometry_stepper.H           Physics/geometry/CD_GeometryStepper.H
# mv Physics/geometry/geometry_stepper.cpp         Physics/geometry/CD_GeometryStepper.cpp

# mv Physics/geometry Physics/Geometry


# for i in `find . -type f \( -iname \*.py -o -iname *GNUmakefile \)`; do
#     sed -i 's/Physics\/geometry/Physics\/Geometry/g' $i
# done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/physics::geometry/Physics::Geometry/g' $i
done
