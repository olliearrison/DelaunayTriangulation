#!/bin/bash

THREADS=(1 2 4 8 10 12)
INPUTS=("inputs/big10000.txt")
MODES=("D" "P" "S" "I")

make clean
make

for input in "${INPUTS[@]}"; do
  echo "===================================="
  echo "Input: $input"
  echo "===================================="

  for mode in "${MODES[@]}"; do
    echo "--- Mode: $mode ---"

    for n in "${THREADS[@]}"; do
      echo "Threads: $n"
      ./triangle -f "$input" -n "$n" -m "$mode" -b 1 -i 10 | grep -E "Parallel mode|Computation time|Valid|Invalid"
    done

    echo
  done
done