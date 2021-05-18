# # Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"advection_diffusion_stepper.H\"/\#include <CD_AdvectionDiffusionStepper.H>/g' $i
#     sed -i 's/\#include \"advection_diffusion_species.H\"/\#include <CD_AdvectionDiffusionSpecies.H>/g' $i
#     sed -i 's/\#include \"advection_diffusion_tagger.H\"/\#include <CD_AdvectionDiffusionTagger.H>/g' $i
# done


# # 
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
#     sed -i 's/advection_diffusion_stepper/AdvectionDiffusionStepper/g' $i
#     sed -i 's/advection_diffusion_species/AdvectionDiffusionSpecies/g' $i
#     sed -i 's/advection_diffusion_tagger/AdvectionDiffusionTagger/g' $i
# done

# for i in `find . -name "GNUmakefile" -type f`;
# do
#     sed -i 's/physics\/advection_diffusion/Physics\/AdvectionDiffusion/g' $i
# done

# for i in `find . -name "*.py" -type f`;
# do
#     sed -i 's/physics\/advection_diffusion/Physics\/AdvectionDiffusion/g' $i
# done

# # Move files
# mv physics/advection_diffusion/advection_diffusion_stepper.H       physics/advection_diffusion/CD_AdvectionDiffusionStepper.H
# mv physics/advection_diffusion/advection_diffusion_stepper.cpp     physics/advection_diffusion/CD_AdvectionDiffusionStepper.cpp
# mv physics/advection_diffusion/advection_diffusion_stepper.options physics/advection_diffusion/CD_AdvectionDiffusionStepper.options
# mv physics/advection_diffusion/advection_diffusion_tagger.H        physics/advection_diffusion/CD_AdvectionDiffusionTagger.H
# mv physics/advection_diffusion/advection_diffusion_tagger.cpp      physics/advection_diffusion/CD_AdvectionDiffusionTagger.cpp
# mv physics/advection_diffusion/advection_diffusion_species.H       physics/advection_diffusion/CD_AdvectionDiffusionSpecies.H
# mv physics/advection_diffusion/advection_diffusion_species.cpp     physics/advection_diffusion/CD_AdvectionDiffusionSpecies.cpp

# mv physics/advection_diffusion Physics/AdvectionDiffusion


for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/advection_diffusion/AdvectionDiffusion/g' $i
done
