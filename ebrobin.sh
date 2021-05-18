# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include <robinconductivityebbcfactory.H>/\#include <CD_RobinConductivityEbBcFactory.H>/g' $i
    sed -i 's/\#include \"robinconductivityebbcfactory.H\"/\#include <CD_RobinConductivityEbBcFactory.H>/g' $i
    
    sed -i 's/\#include <robinconductivityebbc.H>/\#include <CD_RobinConductivityEbBc.H>/g' $i
    sed -i 's/\#include \"robinconductivityebbc.H\"/\#include <CD_RobinConductivityEbBc.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do

    sed -i 's/robinconductivityebbcfactory/RobinConductivityEbBcFactory/g' $i
    sed -i 's/robinconductivityebbc/RobinConductivityEbBc/g' $i
done

# Move files
mv Source/Elliptic/robinconductivityebbc.H     Source/Elliptic/CD_RobinConductivityEbBc.H
mv Source/Elliptic/robinconductivityebbc.cpp   Source/Elliptic/CD_RobinConductivityEbBc.cpp

mv Source/Elliptic/robinconductivityebbcfactory.H     Source/Elliptic/CD_RobinConductivityEbBcFactory.H
mv Source/Elliptic/robinconductivityebbcfactory.cpp   Source/Elliptic/CD_RobinConductivityEbBcFactory.cpp
