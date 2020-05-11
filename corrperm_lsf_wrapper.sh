#!/bin/sh

# set up matlab runtime environment using path stored when module was compiled
export MLR=`cat $2/matlabroot`
#MLR=/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014a

# add MCR paths to library search path
echo setting MCR environment using matlabroot $MLR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MLR/bin/glnxa64:$MLR/sys/java/jre/glnxa64/jre/lib/amd64:$MLR/sys/java/jre/glnxa64/jre/lib/amd64/server:$MLR/sys/java/jre/glnxa64/jre/lib/amd64/jli:$MLR/sys/os/glnxa64:$MLR/runtime/glnxa64

#~old Broad dotkit code
#~  source /broad/software/scripts/useuse
#~  reuse Matlab_2012b_MCR

#run command line
$*
