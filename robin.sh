# # Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include <robinconductivitydomainbcfactory.H>/\#include <CD_RobinConductivityDomainBcFactory.H>/g' $i
#     sed -i 's/\#include \"robinconductivitydomainbcfactory.H\"/\#include <CD_RobinConductivityDomainBcFactory.H>/g' $i
    
#     sed -i 's/\#include <robinconductivitydomainbc.H>/\#include <CD_RobinConductivityDomainBc.H>/g' $i
#     sed -i 's/\#include \"robinconductivitydomainbc.H\"/\#include <CD_RobinConductivityDomainBc.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do

#     sed -i 's/robinconductivitydomainbcfactory/RobinConductivityDomainBcFactory/g' $i
#     sed -i 's/robinconductivitydomainbc/RobinConductivityDomainBc/g' $i
# done

# # Move files
# mv Source/Elliptic/robinconductivitydomainbc.H     Source/Elliptic/CD_RobinConductivityDomainBc.H
# mv Source/Elliptic/robinconductivitydomainbc.cpp   Source/Elliptic/CD_RobinConductivityDomainBc.cpp

# mv Source/Elliptic/robinconductivitydomainbcfactory.H     Source/Elliptic/CD_RobinConductivityDomainBcFactory.H
# mv Source/Elliptic/robinconductivitydomainbcfactory.cpp   Source/Elliptic/CD_RobinConductivityDomainBcFactory.cpp


for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/setCoefficientss/setCoefficients/g' $i
    sed -i 's/m_aCoefficientdata/m_acodata/g' $i
done
