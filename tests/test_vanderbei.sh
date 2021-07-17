#!/usr/bin/env bash
set -euo pipefail

INPFILES=(./test_LPs/input/vanderbei*.txt)
EXPFILES=(./test_LPs/output/vanderbei*.txt)

len=${#INPFILES[@]}

for (( i=0; i<$len; i++))
do
    echo "solving ${INPFILES[$i]}:"
    output=$(../bin/main ${INPFILES[$i]})
    exp_output=$(cat ${EXPFILES[$i]})


    if [[ $(echo "$output" | cut -f1 -d$'\n') = $(echo "$exp_output" | cut -f1 -d$'\n') ]]; then
        echo "passed"
    else
        echo "failed"
        echo "expected:"
        echo "$exp_output"
        echo "got:"
        echo "$output"
    fi

done
