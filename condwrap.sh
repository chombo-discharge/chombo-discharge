# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"conductivitydomainbc_wrapper.H\"/\#include <CD_ConductivityDomainBcWrapper.H>/g' $i
    sed -i 's/\#include \"conductivitydomainbc_wrapper_factory.H\"/\#include <CD_ConductivityDomainBcWrapperFactory.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/conductivitydomainbc_wrapper_factory/ConductivityDomainBcWrapperFactory/g' $i
    sed -i 's/conductivitydomainbc_wrapper/ConductivityDomainBcWrapper/g' $i
    sed -i 's/set_wallbc/setWallBc/g' $i
    sed -i 's/set_potentials/setPotentials/g' $i
    sed -i 's/set_robin_coefs/setRobinCoefficients/g' $i
    sed -i 's/set_coef/setCoefficients/g' $i
    sed -i 's/\.Realm/\.realm/g' $i
    sed -i 's/\.LoadBalancing/\.load_balance/g' $i
    sed -i 's/\.BoxSorting/\.box_sorting/g' $i
done

# Move files
mv Source/elliptic/conductivitydomainbc_wrapper_factory.H   Source/elliptic/CD_ConductivityDomainBcWrapperFactory.H
mv Source/elliptic/conductivitydomainbc_wrapper_factory.cpp Source/elliptic/CD_ConductivityDomainBcWrapperFactory.cpp
mv Source/elliptic/conductivitydomainbc_wrapper.H           Source/elliptic/CD_ConductivityDomainBcWrapper.H
mv Source/elliptic/conductivitydomainbc_wrapper.cpp         Source/elliptic/CD_ConductivityDomainBcWrapper.cpp

# Move elliptic folder to
mv Source/elliptic Source/Elliptic

# Update Py and makefiles
for i in `find . -type f \( -iname \*GNUmakefile -o -iname \*.py \)`; do
    sed -i 's/elliptic/Elliptic/g' $i
done
