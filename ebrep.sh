for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/EB_representation::implicit_function/EB_representation::ImplicitFunction/g' $i
    sed -i 's/EB_representation::discrete/EB_representation::Discrete/g' $i
    sed -i 's/EB_representation::voxel/EB_representation::Voxel/g' $i
    sed -i 's/EB_representation::level_set/EB_representation::Levelset/g' $i
    sed -i 's/EB_representation/EbRepresentation/g' $i
done

