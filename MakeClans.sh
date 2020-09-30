#!/usr/bin/bash

rm -f -r clans/
javac -d ./ src/clans/*.java \
            src/clans/model/*.java \
            src/clans/model/proteins/*.java \
            src/clans/model/microarray/*.java \
            src/clans/misc/*.java \
            src/clans/io/*.java \
            src/clans/headless/*.java \
            src/clans/gui/*.java \
            src/clans/algorithms/*.java \
            src/clans/algorithms/fruchtermanreingold/*.java

jar cfe clans.jar clans.Main clans/
