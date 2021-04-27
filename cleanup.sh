base=field_solver
BASE=FIELD_SOLVER
implem=field_solver_multigrid
IMPLEM=FIELD_SOLVER_MULTIGRID
factory=field_solver_factory
factoryI=field_solver_factoryI

# Move poisson_solver.* field_solver.*
# mv src/field_solver/poisson_solver.cpp src/field_solver/$base.cpp
# mv src/field_solver/poisson_solver.H src/field_solver/$base.H
# mv src/field_solver/poisson_solver.options src/field_solver/$base.options

# # Same for implementation class
# mv src/field_solver/poisson_multifluid_gmg.cpp src/field_solver/$implem.cpp
# mv src/field_solver/poisson_multifluid_gmg.H src/field_solver/$implem.H
# mv src/field_solver/poisson_multifluid_gmg.options src/field_solver/$implem.options

# # Factory class
# mv src/field_solver/poisson_factory.H src/field_solver/$factory.H
# mv src/field_solver/poisson_factoryI.H src/field_solver/$factoryI.H

# Move poisson_stepper to field_stepper
mv physics/poisson/poisson_stepper.H physics/poisson/field_stepper.H
mv physics/poisson/poisson_stepperI.H physics/poisson/field_stepperI.H
mv physics/poisson/poisson_stepper.options physics/poisson/field_stepper.options
mv physics/poisson physics/field

# Move regression tests
#mv regression/poisson regression/field
#mv regression/poisson.ini regression/field.ini


# Find and replaces field_solver with field_solver stuff in all .H files
for i in `find . -name "*.H" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
    sed -i 's/POISSON_SOLVER/FIELD_SOLVER/g' $i
    sed -i 's/poisson_stepper/field_stepper/g' $i
    sed -i 's/POISSON_STEPPER/FIELD_STEPPER/g' $i
done
for i in `find . -name "*.cpp" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
    sed -i 's/POISSON_SOLVER/FIELD_SOLVER/g' $i
    sed -i 's/poisson_stepper/field_stepper/g' $i
    sed -i 's/POISSON_STEPPER/FIELD_STEPPER/g' $i    
done
for i in `find . -name "*.options" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
    sed -i 's/POISSON_SOLVER/FIELD_SOLVER/g' $i
    sed -i 's/poisson_stepper/field_stepper/g' $i
    sed -i 's/POISSON_STEPPER/FIELD_STEPPER/g' $i    
done
for i in `find . -name "*.inputs" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
    sed -i 's/POISSON_SOLVER/FIELD_SOLVER/g' $i
    sed -i 's/poisson_stepper/field_stepper/g' $i
    sed -i 's/POISSON_STEPPER/FIELD_STEPPER/g' $i    
done
for i in `find . -name "*.py" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
    sed -i 's/POISSON_SOLVER/FIELD_SOLVER/g' $i
    sed -i 's/poisson_stepper/field_stepper/g' $i
    sed -i 's/POISSON_STEPPER/FIELD_STEPPER/g' $i    
done
for i in `find . -name "GNUmakefile" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
    sed -i 's/POISSON_SOLVER/FIELD_SOLVER/g' $i
    sed -i 's/poisson_stepper/field_stepper/g' $i
    sed -i 's/POISSON_STEPPER/FIELD_STEPPER/g' $i
    sed -i 's/poisson/field/g' $i        
done

# Find and replaces poisson_multifluid_gmg with field_solver stuff in all .H files
for i in `find . -name "*.H" -type f`;
do
    sed -i "s/poisson_multifluid_gmg/$implem/g" $i
    sed -i "s/POISSON_MULTIFLUID_GMG/$IMPLEM/g" $i    
done
for i in `find . -name "*.cpp" -type f`;
do
    sed -i "s/poisson_multifluid_gmg/$implem/g" $i
    sed -i "s/POISSON_MULTIFLUID_GMG/$IMPLEM/g" $i    
done
for i in `find . -name "*.options" -type f`;
do
    sed -i "s/poisson_multifluid_gmg/$implem/g" $i
    sed -i "s/POISSON_MULTIFLUID_GMG/$IMPLEM/g" $i        
done
for i in `find . -name "*.inputs" -type f`;
do
    sed -i "s/poisson_multifluid_gmg/$implem/g" $i
    sed -i "s/POISSON_MULTIFLUID_GMG/$IMPLEM/g" $i        
done
for i in `find . -name "*.py" -type f`;
do
    sed -i "s/poisson_multifluid_gmg/$implem/g" $i
    sed -i "s/POISSON_MULTIFLUID_GMG/$IMPLEM/g" $i        
done

# Factory rename
# for i in `find . -name "*.H" -type f`;
# do
#     sed -i "s/poisson_factory/$factory/g" $i
#     sed -i "s/poisson_factoryI/$factoryI/g" $i
# done
# for i in `find . -name "*.cpp" -type f`;
# do
#     sed -i "s/poisson_factory/$factory/g" $i
#     sed -i "s/poisson_factoryI/$factoryI/g" $i    
# done

# Fix some typos in input and option files. 
# for i in `find . -name "*.options" -type f`;
# do
#     sed -i "s/GMG_GMG/GMG/g" $i
# done
# for i in `find . -name "*.inputs" -type f`;
# do
#     sed -i "s/GMG_GMG/GMG/g" $i
# done    

# Replace regtest names
for i in `find . -name "*.ini" -type f`;
do
    sed -i 's/poisson/field/g' $i
done
