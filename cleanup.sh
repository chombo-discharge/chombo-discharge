base=field_solver
BASE=field_solver
implem=field_solver_multigrid
IMPLEM=FIELD_SOLVER_MULTIGRID
factory=field_solver_factory
factoryI=field_solver_factoryI

echo $newName

# Move poisson_solver.* field_solver.*
mv src/field_solver/poisson_solver.cpp src/field_solver/$base.cpp
mv src/field_solver/poisson_solver.H src/field_solver/$base.H
mv src/field_solver/poisson_solver.options src/field_solver/$base.options

# Same for implementation class
mv src/field_solver/poisson_multifluid_gmg.cpp src/field_solver/$implem.cpp
mv src/field_solver/poisson_multifluid_gmg.H src/field_solver/$implem.H
mv src/field_solver/poisson_multifluid_gmg.options src/field_solver/$implem.options

# Factory class
mv src/field_solver/poisson_factory.H src/field_solver/$factory.H
mv src/field_solver/poisson_factoryI.H src/field_solver/$factoryI.H

# Find and replaces field_solver with field_solver stuff in all .H files
for i in `find . -name "*.H" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
done
for i in `find . -name "*.cpp" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
done
for i in `find . -name "*.options" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
done
for i in `find . -name "*.inputs" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
done
for i in `find . -name "*.py" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
done
for i in `find . -name "GNUmakefile" -type f`;
do
    sed -i 's/poisson_solver/field_solver/g' $i
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

#
for i in `find . -name "*.H" -type f`;
do
    sed -i "s/poisson_factory/$factory/g" $i
done
for i in `find . -name "*.cpp" -type f`;
do
    sed -i "s/poisson_factory/$factory/g" $i
done

# Fix some typos in input and option files. 
for i in `find . -name "*.options" -type f`;
do
    sed -i "s/GMG_GMG/GMG/g" $i
done
for i in `find . -name "*.inputs" -type f`;
do
    sed -i "s/GMG_GMG/GMG/g" $i
done    
