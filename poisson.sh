base=field_solver
BASE=FIELD_SOLVER
implem=field_solver_multigrid
IMPLEM=FIELD_SOLVER_MULTIGRID
factory=field_solver_factory
factoryI=field_solver_factoryI

# Factory rename
for i in `find . -name "*.H" -type f`;
do
    sed -i "s/poisson_factory/$factory/g" $i
    sed -i "s/poisson_factoryI/$factoryI/g" $i
done
for i in `find . -name "*.cpp" -type f`;
do
    sed -i "s/poisson_factory/$factory/g" $i
    sed -i "s/poisson_factoryI/$factoryI/g" $i    
done
for i in `find . -name "*.py" -type f`;
do
    sed -i "s/poisson_factory/$factory/g" $i
    sed -i "s/poisson_factoryI/$factoryI/g" $i    
done

