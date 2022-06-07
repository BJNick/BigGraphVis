# A shell script to run ./graph_viewer script over and over with the same parameters
# Prints complete at the end
#
# Usage:
#   ./repeat.sh <script> <parameters>
#
# Example:
#   ./repeat.sh lucene
#

for i in {0..4}
do
    ./graph_viewer $1 $2 $3 $4 $5 $6 $7 $8 $9 rand$i
    echo "█ FINISHED ITERATION RAND $i █"
done

. ~/gpu2/bigtext.sh