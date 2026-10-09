#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
HOSTFILE="MPI_Local_Host"
LCN=4
STATE=2
DATAFILE="NX360_NP128"
INPUT=${1:-Data/${LCN}/${STATE}/${DATAFILE}}
OUTPUT=${2:-${INPUT}_output}
mpirun --hostfile $HOSTFILE ./build/main $LCN $INPUT 1E-2 1E-5 2 0.71 1 > $OUTPUT &
