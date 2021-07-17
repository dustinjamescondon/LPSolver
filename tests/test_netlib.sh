#!/usr/bin/env bash
set -euo pipefail

INPFILES=(./tests/test_LPs/input/netlib*.txt)
EXPFILES=(./tests/test_LPs/output/netlib*.txt)

len=${#INPFILES[@]}

for (( i=0; i<$len; i++))
do
    output=$(./bin/main ${INPFILES[$i]})
    exp_output=$(cat ${EXPFILES[$i]})

    passed=false
    if [[ $(echo "$output" | cut -f1 -d$'\n') = $(echo "$exp_output" | cut -f1 -d$'\n') ]]; then
        if [[ $(echo "$output" | cut -f2 -d$'\n') = $(echo "$exp_output" | cut -f2 -d$'\n') ]]; then
            passed=true
        fi
    fi

    if [ "$passed" = true ]; then
        echo "passed"
    else
        echo "--------------------------------------------------"
        echo "failed"
        echo ".................................................."
        echo "solving ${INPFILES[$i]}:"
        echo "expected:"
        echo "$exp_output"
        echo "got:"
        echo "$output"
        echo "--------------------------------------------------"
    fi
done
