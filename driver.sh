# Delete flags from input files
for i in `find . -name "*.inputs" -type f`; do
    sed -i '/driver.recursive_regrid/d' $i
    sed -i '/driver.write_ebis/d' $i
    sed -i '/driver.read_ebis/d' $i
    sed -i '/driver.grow_tags/d' $i
done

# Delete flags from option files
for i in `find . -name "*.options" -type f`; do
    sed -i '/driver.recursive_regrid/d' $i
    sed -i '/driver.write_ebis/d' $i
    sed -i '/driver.read_ebis/d' $i
    sed -i '/driver.grow_tags/d' $i
done

# Find and replaces setup_and_run() stuff in all main.cpp files
for i in `find . -name "main.cpp" -type f`;
do
    sed -i 's/char\* input_file = argv\[1\]\;/const std::string input_file = argv[1];/g' $i
    sed -i 's/ParmParse pp(argc-2, argv+2, NULL, input_file);/ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());/g' $i
    sed -i 's/engine->setup_and_run()/engine->setup_and_run(input_file)/g' $i
done

# Do the same in all Python files
for i in `find . -name "*.py" -type f`;
do
    sed -i 's/char\* input_file = argv\[1\]\;/const std::string input_file = argv[1];/g' $i
    sed -i 's/ParmParse pp(argc-2, argv+2, NULL, input_file);/ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());/g' $i
    sed -i 's/engine->setup_and_run()/engine->setup_and_run(input_file)/g' $i
done

# Insert driver.write_restart_file into all input scripts
for i in `find . -name "*.inputs" -type f`;
do
    lineNumber=`grep -n "driver.write_regrid_files"  $i | cut -f1 -d:`
    if [ -z "$lineNumber" ]
    then
	echo $i
    else
	sed -i "$lineNumber a driver.write_restart_files             = false         # Write restart files or not" $i
    fi
done

# Find all write_main.py and delete them
for i in `find . -name "write_main.py" -type f`;
do
    rm $i
done
