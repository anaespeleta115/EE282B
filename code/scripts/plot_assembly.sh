#!/usr/bin/env bash
set -euo pipefail

OUTDIR="./out"
PLOTCDF2="/data/homezvol2/aespelet/bin/plotCDF2"

mkdir -p "$OUTDIR"

# Build plotCDF2-formatted input for your assembly
bioawk -c fastx '{ print length($seq) }' "$OUTDIR/hifi_fly_assembly.bp.p_ctg.fa" \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nMy_Assembly\t0" } { print "My_Assembly\t" $1 }' \
> "$OUTDIR/my_assembly_cdf2.txt"

# Build plotCDF2-formatted input for FlyBase scaffolds
bioawk -c fastx '{ print length($seq) }' "$OUTDIR/flybase_scaffolds.fa" \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nFB_Scaff\t0" } { print "FB_Scaff\t" $1 }' \
> "$OUTDIR/flybase_scaffolds_cdf2.txt"

# Build plotCDF2-formatted input for FlyBase contigs
bioawk -c fastx '{ print length($seq) }' "$OUTDIR/flybase_contigs.fa" \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nFB_Ctg\t0" } { print "FB_Ctg\t" $1 }' \
> "$OUTDIR/flybase_contigs_cdf2.txt"

# Run plotCDF2
"$PLOTCDF2" \
"$OUTDIR/flybase_scaffolds_cdf2.txt" \
"$OUTDIR/flybase_contigs_cdf2.txt" \
"$OUTDIR/my_assembly_cdf2.txt" \
"$OUTDIR/contiguity_plot2.png"

echo "Done: $OUTDIR/contiguity_plot.png"
