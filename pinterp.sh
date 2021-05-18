#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"EBParticleInterp.H\"/\#include <CD_EbParticleInterp.H>/g' $i
    sed -i 's/\#include <EBParticleInterp.H>/\#include <CD_EbParticleInterp.H>/g' $i
    sed -i 's/\#include \"EBParticleInterpI.H\"/\#include <CD_EbParticleInterpImplem.H>/g' $i
    sed -i 's/\#include <EBParticleInterpI.H>/\#include <CD_EbParticleInterpImplem.H>/g' $i
    sed -i 's/\#include <EBParticleInterpF_F.H>/\#include <CD_EbParticleInterpF_F.H>/g' $i
    sed -i 's/\#include \"EBParticleInterpF_F.H\"/\#include <CD_EbParticleInterpF_F.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/EBParticleInterp/EbParticleInterp/g' $i
done

# Move files
mv Source/Particle/EBParticleInterp.H     Source/Particle/CD_EbParticleInterp.H
mv Source/Particle/EBParticleInterpI.H    Source/Particle/CD_EbParticleInterpImplem.H
mv Source/Particle/EBParticleInterp.cpp   Source/Particle/CD_EbParticleInterp.cpp
mv Source/Particle/EBParticleInterpF_F.H  Source/Particle/CD_EbParticleInterpF_F.H
mv Source/Particle/EBParticleInterpF.ChF  Source/Particle/CD_EbParticleInterpF.ChF


