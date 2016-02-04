#!/bin/bash
export PATH=/kb/runtime/bin:/mnt/kb/dev_container/bin:/kb/runtime/glassfish3/bin:/kb/runtime//java/bin:/kb/runtime/thrift/bin:/kb/runtime/bin:/mnt/kb/dev_container/bin:/kb/runtime/glassfish3/bin:/kb/runtime//java/bin:/kb/runtime/thrift/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/kb/deployment/bin
export KB_RUNTIME=/kb/runtime
export PYTHONPATH="/mnt/kb/dev_container/modules/coexpression/lib"
python /mnt/kb/dev_container/modules/coexpression/ltest/script_test/basic_test.py $1 $2 $3
