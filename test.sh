#!/bin/bash
for f in $(find ./regression/ -type f -name '*.inputs')
do
    sed -i '21 i amr_mesh.lsf_ghost       = 3           # Number of ghost cells when writing level-set to grid' $f
done
