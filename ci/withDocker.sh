#!/bin/bash
CI_DIR=$(cd "$(dirname "$0")"; pwd)
PBBAMR_DIR=$(cd "$(dirname "$0")"/..; pwd)
docker build -t pbbamr $CI_DIR
docker run -h docker-host -v $PBBAMR_DIR:/pbbamr -w /pbbamr --rm=true  -i -t pbbamr $*
