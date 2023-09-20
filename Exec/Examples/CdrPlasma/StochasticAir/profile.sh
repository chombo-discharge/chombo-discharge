for PATTERN in 'McPhoto::advancePhotonsInstantaneous::amr_loop' \
		   'McPhoto::advancePhotonsInstantaneous::remap_bulk' \
		   'McPhoto::advancePhotonsInstantaneous::mpi_waitall' \
	       ; do

    for i in `find . -name "time.table.*"`; do	
	grep -rn -e $PATTERN $i | head -1
    done

    echo "";
done
