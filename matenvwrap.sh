#!/bin/sh
# set up matlab runtime environment using path stored when module was compiled
export MLR=`cat $2/matlabroot`
# add MCR paths to library search path
echo setting MCR environment using matlabroot $MLR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MLR/bin/glnxa64:$MLR/sys/java/jre/glnxa64/jre/lib/amd64:$MLR/sys/java/jre/glnxa64/jre/lib/amd64/server:$MLR/sys/java/jre/glnxa64/jre/lib/amd64/jli:$MLR/sys/os/glnxa64:$MLR/runtime/glnxa64
#run command line
$*


