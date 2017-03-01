#!/bin/bash
# ---
# This script is called by the bamboo build job after a successful
# build on the master branch.  It deploys only pbbamr.
#
# Note that the build of the underlying R stack---all the third
# party R libraries---is handled elsewhere:
# - http://jenkins:8080/view/ITG/job/InternalTools_R_Build_CentOS7/
# - http://jenkins:8080/view/ITG/job/InternalTools_R_Build_Ubuntu14/
# ---


# Expects an argument centos7 or u1404 to

dist=$1

source /mnt/software/Modules/current/init/bash
module use /mnt/software/modulefiles
module purge
module load R/3.2.2-internal

case $dist in
    centos7)
        export R_LIBS=/pbi/dept/secondary/builds/R/current/lib/R/library/centos7
        echo "Deploying for CentOS7"
        ;;
    u1404)
        export R_LIBS=/pbi/dept/secondary/builds/R/current/lib/R/library/u1404
        echo "Deploying for Ubuntu 14.04"
        ;;
    *)
        echo "Invalid or missing argument: must specify centos7 or u1404 as single argument"
        exit 1
        ;;
esac

# Remove any object files
find . -name "*.o" -o -name "*.so"  -o -name "*.a" | xargs rm


# Install pbbamr
R CMD INSTALL .
