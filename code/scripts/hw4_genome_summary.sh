#!/usr/bin/env bash
set -euo pipefail

GENOME="../../data/dmel-all-chromosome-r6.66.fasta.gz"
OUTDIR="./out"
mkdir -p "$OUTDIR"

gzip -dc "$GENOME" | bioawk -c fastx '
BEGIN {
    small_bases = small_Ns = small_seqs = 0
    large_bases = large_Ns = large_seqs = 0
}
{
    len = length($seq)

    seq_copy = $seq
    ncount = gsub(/[Nn]/, "", seq_copy)

    if (len <= 100000) {
        small_bases += len
        small_Ns += ncount
        small_seqs += 1
    } else {
        large_bases += len
        large_Ns += ncount
        large_seqs += 1
    }
}
END {
    print "partition\tmetric\tvalue"
    print "<=100kb\ttotal_nucleotides\t" small_bases
    print "<=100kb\ttotal_Ns\t" small_Ns
    print "<=100kb\ttotal_sequences\t" small_seqs
    print ">100kb\ttotal_nucleotides\t" large_bases
    print ">100kb\ttotal_Ns\t" large_Ns
    print ">100kb\ttotal_sequences\t" large_seqs
}
' > "$OUTDIR/genome_summary.tsv"
