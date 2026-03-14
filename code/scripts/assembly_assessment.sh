#!/usr/bin/env bash

set -euo pipefail

GFA="./out/hifi_fly_assembly.bp.p_ctg.gfa"
FA="./out/hifi_fly_assembly.bp.p_ctg.fa"



awk '/^S/{print ">"$2; print $3}' "$GFA" > "$FA"



bioawk -c fastx '{print length($seq)}' "$FA" \
| sort -nr \
> ./out/contig_lengths_sorted.txt



TOTAL=$(awk '{sum+=$1} END{print sum}' ./out/contig_lengths_sorted.txt)

echo "Total assembly size: $TOTAL"


N50=$(awk -v half=$((TOTAL/2)) '
{
    cum += $1
    if (cum >= half) {
        print $1
        exit
    }
}' ./out/contig_lengths_sorted.txt)

echo "N50: $N50"
