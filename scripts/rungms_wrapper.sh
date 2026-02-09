#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 1 ]; then
  echo "usage: $0 <job.inp>" >&2
  exit 2
fi
input="$1"
job="${input%.inp}"
job="$(basename "$job")"
# VERNO=00, NCPUS=1, PPN=1, LOGN=0
exec /home/fkulcs/molcalc/scripts/rungms_csh "$job" 00 1 1 0
