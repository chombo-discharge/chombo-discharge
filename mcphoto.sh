# #Inclusion guards
# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
#     sed -i 's/\#include \"mc_photo.H\"/\#include <CD_McPhoto.H>/g' $i
#     sed -i 's/\#include <mc_photo.H>/\#include <CD_McPhoto.H>/g' $i
# done

# for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
#     sed -i 's/mc_photo/McPhoto/g' $i
#     sed -i 's/is_instantaneous/isInstantaneous/g' $i
#     sed -i 's/generate_Photons/generatePhotons/g' $i
#     sed -i 's/advance_Photons_stationary/advancePhotonsStationary/g' $i
#     sed -i 's/advance_Photons_transient/advancePhotonsTransient/g' $i
#     sed -i 's/deposit_Photons/depositPhotons/g' $i
#     sed -i 's/sort_Photons_by_cell/sortPhotonsByCell/g' $i
#     sed -i 's/sort_Photons_by_patch/sortPhotonsByPatch/g' $i
#     sed -i 's/sort_bulk_Photons_by_cell/sortBulkPhotonsByCell/g' $i
#     sed -i 's/sort_bulk_Photons_by_patch/sortBulkPhotonsByPatch/g' $i
#     sed -i 's/sort_source_Photons_by_cell/sortSourcePhotonsByCell/g' $i
#     sed -i 's/sort_source_Photons_by_patch/sortSourcePhotonsByPatch/g' $i
#     sed -i 's/sort_eb_Photons_by_cell/sortEbPhotonsByCell/g' $i
#     sed -i 's/sort_eb_Photons_by_patch/sortEbPhotonsByPatch/g' $i
#     sed -i 's/sort_domain_Photons_by_cell/sortDomainPhotonsByCell/g' $i
#     sed -i 's/sort_domain_Photons_by_patch/sortDomainPhotonsByPatch/g' $i
#     sed -i 's/getPVR_buffer/getPVRBuffer/g' $i
#     sed -i 's/get_halo_buffer/getHaloBuffer/g' $i
#     sed -i 's/set_pvr_buffer/setPVRBuffer/g' $i
#     sed -i 's/set_halo_buffer/setHalobuffer/g' $i
#     sed -i 's/count_Photons/countPhotons/g' $i
#     sed -i 's/count_outcast/countOutcast/g' $i
#     sed -i 's/get_Photons/getPhotons/g' $i
#     sed -i 's/get_bulk_Photons/getBulkPhotons/g' $i
#     sed -i 's/get_eb_Photons/getEbPhotons/g' $i
#     sed -i 's/get_domain_Photons/getDomainPhotons/g' $i
#     sed -i 's/getSource_Photons/getSourcePhotons/g' $i
#     sed -i 's/Photon_generation/PhotonGeneration/g' $i
#     sed -i 's/Source_type/SourceType/g' $i
#     sed -i 's/draw_photons/drawPhotons/g' $i
#     sed -i 's/domainbc_map/domainBcMap/g' $i
#     sed -i 's/random_fieldSolver/randomPoisson/g' $i
#     sed -i 's/random_exponential/randomExponential/g' $i
#     sed -i 's/random_direction/randomDirection/g' $i
#     sed -i 's/deposit_kappaConservative/depositKappaConservative/g' $i
#     sed -i 's/deposit_nonConservative/depositNonConservative/g' $i
#     sed -i 's/deposit_hybrid/depositHybrid/g' $i
#     sed -i 's/level_redistribution/levelRedist/g' $i
#     sed -i 's/parse_rng/parseRng/g' $i
#     sed -i 's/parse_pseudoPhotons/parsePseudoPhotons/g' $i
#     sed -i 's/parse_photogen/parsePhotoGeneration/g' $i
#     sed -i 's/parse_instantaneous/parseInstantaneous/g' $i
#     sed -i 's/parse_source_type/parseSourceType/g' $i
#     sed -i 's/parse_deposition/parseDeposition/g' $i
#     sed -i 's/parse_bisect_step/parseBisectStep/g' $i
#     sed -i 's/parse_pvr_buffer/parsePvrBuffer/g' $i
# done

# # Move files
# mv Source/RadiativeTransfer/mc_photo.H       Source/RadiativeTransfer/CD_McPhoto.H
# mv Source/RadiativeTransfer/mc_photo.cpp     Source/RadiativeTransfer/CD_McPhoto.cpp
# mv Source/RadiativeTransfer/mc_photo.options Source/RadiativeTransfer/CD_McPhoto.options

for i in `find . -type f \( -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/.max_Photons/.max_photons/g' $i
    sed -i 's/.PhotonGeneration/.photon_generation/g' $i
done
