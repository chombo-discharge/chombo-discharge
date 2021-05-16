# Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"cell_tagger.H\"/\#include \"CD_CellTagger.H\"/g' $i
    sed -i 's/\#include <cell_tagger.H>/\#include <CD_CellTagger.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cell_tagger/CellTagger/g' $i

done

# Move files
mv src/cell_tagger/cell_tagger.cpp     src/cell_tagger/CD_CellTagger.cpp
mv src/cell_tagger/cell_tagger.H       src/cell_tagger/CD_CellTagger.H
mv src/cell_tagger/cell_tagger.options src/cell_tagger/CD_CellTagger.options

# Move folder
mv src/cell_tagger src/CellTagger

# 
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options \)`; do
    sed -i 's/cell_tagger/CellTagger/g' $i
    sed -i 's/parse_boxes/parseBoxes/g' $i
    sed -i 's/parse_buffer/parseBuffer/g' $i
    sed -i 's/inside_tag_box/insideTagBox/g' $i
    sed -i 's/m_tagboxes/m_tagBoxes/g' $i
    sed -i 's/get_buffer/getBuffer/g' $i
done

# Update makefiles
for i in `find . -name "*GNUmakefile" -type f`; do
    sed -i 's/cell_tagger/CellTagger/g' $i
    sed -i 's/CellTaggers/cell_taggers/g' $i
done
