#!/usr/bin/env bash
set -euo pipefail

GENOME="../../data/dmel-all-chromosome-r6.66.fasta.gz"
OUTDIR="./out"
mkdir -p "$OUTDIR"

gzip -dc "$GENOME" | bioawk -c fastx '
BEGIN {
    print "seq_name\tseq_length\tgc_percent\tpartition"
}
{
    len = length($seq)
    gc_pct = gc($seq)

    part = (len <= 100000 ? "<=100kb" : ">100kb")

    print $name "\t" len "\t" gc_pct "\t" part
}
' > "$OUTDIR/sequence_stats.tsv"
