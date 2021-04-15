for i in `find . -name "*.inputs" -type f`; do
    echo $i

    # Remove from amr_mesh
    sed -i '/amr_mesh.refine_all_lvl/d' $i
done

for i in `find . -name "*.options" -type f`; do
    echo $i

    # Remove from amr_mesh
    sed -i '/amr_mesh.refine_all_lvl/d' $i
done

