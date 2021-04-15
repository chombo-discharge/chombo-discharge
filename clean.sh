for i in `find . -name "*.options" -type f`; do
    echo $i
     sed -i '/driver.recursive_regrid/d' $i
     sed -i '/driver.write_ebis/d' $i
     sed -i '/driver.read_ebis/d' $i
     sed -i '/driver.grow_tags/d' $i
done

