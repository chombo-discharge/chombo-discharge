# This script find the line number partially matching an expression in a group of files. It then
# append a new expression directly after the match.
#
# Pattern to match and expression to insert
PATTERN='StreamerInception'
REPLACEMENT='DischargeInception'

for i in `find . -type f -name "*.rst" -not -path "*./Submodules/*"`; do
    sed -i "s/$PATTERN/$REPLACEMENT/g" $i
done
