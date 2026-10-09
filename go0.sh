#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
HOSTFILE="MPI_Local_Host"
LCN=6
STATE=2
DATAFILE="NX80_NP6144"
INPUT=${1:-Data/${LCN}/${STATE}/${DATAFILE}}
OUTPUT=${2:-${INPUT}_output}
mpirun --hostfile $HOSTFILE ./build/main $LCN $INPUT 1E-2 1E-5 2 0.25 1 > $OUTPUT &
