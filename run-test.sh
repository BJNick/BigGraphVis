# This script compiles the BigGraphVis project and runs a simple test.
. command.sh
cd builds/linux
make graph_viewer && ./graph_viewer gpu 500 1 sg 80 1 approximate ../../../net/web-BerkStan.txt  ../../../output png 1024 1024 11 5 6500 
cd ../../
