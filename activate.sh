#!/bin/bash
source ../BaseStack/bin/setup_run.sh
PYTHONPATH=`pwd`/src:${PYTHONPATH}:../iBuildSB/src
export PYTHONPATH
source ibuild/bin/activate
