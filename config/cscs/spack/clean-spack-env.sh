#!/bin/bash
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)

echo "Removing all spack environement in $SCRIPT_DIR"
rm  -f "$SCRIPT_DIR"/v*/*/spack.lock
rm -rf "$SCRIPT_DIR"/v*/*/.spack-env
