#!/usr/bin/env bash
set -euo pipefail

INPFILES=(./test_LPs/input/vanderbei*.txt)
EXPFILES=(./test_LPs/output/vanderbei*.txt)

len=${#INPFILES[@]}

for (( i=0; i<$len; i++))
do
    echo "solving ${INPFILES[$i]}:"
    output=$(../bin/main ${INPFILES[$i]})
    echo "$output"

    exp_output=$(cat ${EXPFILES[$i]})
    echo "$exp_output"
done
