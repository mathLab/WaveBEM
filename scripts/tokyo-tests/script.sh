#!/bin/sh

F="0250 0289 0354 0267 0316 0408"

test -f launch.sh && rm launch.sh
for a in $F; do 
    echo "Processing Freude n. $a"
    mydir=tokio_FR$a
    test -d tokio_FR$a || mkdir $mydir
    cd $mydir
    sed "s/XXXX/$a/" ../base.prm > waveBem.prm
    test -L waveBem || ln -s ../../waveBem .
    cp ../base.pbs tokio_FR$a.pbs 
    cd ..
    echo "cd $mydir; qsub tokio_FR$a.pbs; cd .." >> launch.sh
done
