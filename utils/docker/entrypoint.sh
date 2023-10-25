#!/usr/bin/bash

set -x
set -euo pipefail

first=${1-}

if [ "$first" == exec ]; then
  shift
  exec "$@"
else
  exec \
    varfish-server-worker "$@"
fi

exit $?
