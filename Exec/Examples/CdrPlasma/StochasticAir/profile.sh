
for PATTERN in 'ParticleData::remapOutcast' \
		   'ParticleData::remapOutcast::entry_barrier' \
		   'ParticleData::remapOutcast::grid_loop1' \
		   'ParticleData::remapOutcast::rotate' \
		   'ParticleData::remapOutcast::build_map' \
		   'ParticleData::remapOutcast::collect_map' \
		   'ParticleData::remapOutcast::barrier' \
		   'ParticleData::remapOutcast::mpi_scatter' \
		   'ParticleData::remapOutcast::catenate' \
		   'McPhoto::depositPhotons'\
	       ; do

    for i in `find . -name "time.table.*"`; do	
	grep -rn -e $PATTERN $i | head -1
    done

    echo "";
done
