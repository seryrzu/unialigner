#!/usr/bin/env bash

# Usage ./bed_getfasta.sh <input_fasta> <bed> <output_fasta>
bedtools getfasta -fi $1 -bed $2 | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | tr -d '\n' | cat <(head -n 1 $1) - | fold -w 60 > $3
