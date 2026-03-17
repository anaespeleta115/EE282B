---
title: "hw4"
output: md_document
date: "2026-03-18"
author: "Ana Espeleta"
---


## Homework 4:

This document serves as the master report for Homework 4.  
All commands were executed inside my EE282B GitHub repository on the homework4 branch:

### Part 1: Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:

Inside the /code/scripts/ directory, run the following script, which I've named: hw4_genome_summary.sh

```{bash, eval=FALSE}
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

```
From here, I got the following output, which is stored inside /code/scripts/out/genome/summary.tsv:

partition	    metric	      value
<=100kb	  total_nucleotides	6178042
<=100kb	  total_Ns	        662593
<=100kb	  total_sequences	  1863
>100kb	  total_nucleotides	137547960
>100kb	  total_Ns	        490385
>100kb	  total_sequences	  7

### Part 2: Plot the length and GC content distributions

#### Calculate the length and GC content using bioawk:

```{bash, eval=FALSE}
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
```

#### Import sequence stats file into R and plot histograms:


```{R, eval=FALSE}
library(ggplot2)
library(readr)

df <- read_tsv("out/sequence_stats.tsv")

ggplot(df, aes(x = seq_length)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  facet_wrap(~ partition, scales = "free_y") +
  labs(
    title = "Sequence length distribution",
    x = "Sequence length (log scale)",
    y = "Number of sequences"
  ) +
  theme_bw()



ggplot(df, aes(x = gc_percent)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ partition, scales = "free_y") +
  labs(
    title = "Sequence GC% distribution",
    x = "GC percentage",
    y = "Number of sequences"
  ) +
  theme_bw()

```

#### Use plotCDF script to plot cumulative sequence size:

```{bash, eval=FALSE}
awk -F'\t' 'NR>1 && $4=="<=100kb" {print $2}' out/sequence_stats.tsv \
| sort -nr \
> out/lengths_le_100kb.txt

awk -F'\t' 'NR>1 && $4==">100kb" {print $2}' out/sequence_stats.tsv \
| sort -nr \
> out/lengths_gt_100kb.txt
```
plotCDF:

```{bash, eval=FALSE}
/data/homezvol2/aespelet/bin/plotCDF \
out/lengths_le_100kb.txt \
out/cdf_le_100kb.png

/data/homezvol2/aespelet/bin/plotCDF \
out/lengths_gt_100kb.txt \
out/cdf_gt_100kb.png
```

These plots are all stored inside /code/scripts/out/, under the following names:

1. seq_gc_distribution.png (file has two plots)
2. seq_length_distribution_plot.png (file has 2 plots)
3. cdf_gt_100kb.png
4. cdf_le_100kb.png



### Part 3: Assemble a genome using Pacbio HiFi reads

After downloading the reads, I made a script called run_hifiasm.sh

```{bash, eval=FALSE}
!/bin/bash
#SBATCH --account=ecoevo282_class
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err

set -euo pipefail

# Initialize conda
source /dfs6b/pub/aespelet/miniforge3/etc/profile.d/conda.sh

conda activate ee282

# Run hifiasm
hifiasm \
    -o hifi_fly_assembly \
    -t ${SLURM_CPUS_PER_TASK} \
    /pub/jje/ee282/ISO_HiFi_Shukla2025.fasta.gz
```


I ran it using the sbatch command:

```{bash, eval=FALSE}
sbatch ./run_hifiasm
```

### Part 4: Assembly assessment

#### Calculate N50s:

First, I converted the .gfa file to a fasta file:

```{bash, eval=FALSE}
awk '/^S/{print ">"$2; print $3}' hifi_fly_assembly.bp.p_ctg.gfa > hifi_fly_assembly.bp.p_ctg.fa
```
Next, I made a script called assembly_assessment.sh, which gets the contig lengths, sorts lengths from smallest to largest, and computes the N50. 

After running the assessment script, I got that my N50 was 21.7 Mb, while the reference's was 21.5 Mb.

#### Plot the assembly assessment:

To plot the assessment, I made a new script, plot_assembly.sh, to extract sequence lengths from the three assemblies (my assembly and FlyBase's contig and scaffold assemblies), format them for plotCDF2, and generate the contiguity plot.

```{bash, eval=FALSE}
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="./out"
PLOTCDF2="/data/homezvol2/aespelet/bin/plotCDF2"

mkdir -p "$OUTDIR"

# Build plotCDF2-formatted input for my assembly
bioawk -c fastx '{ print length($seq) }' "$OUTDIR/hifi_fly_assembly.bp.p_ctg.fa>
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nMy_Assembly\t0" } { print "My_Assembly\>
> "$OUTDIR/my_assembly_cdf2.txt"

# Build plotCDF2-formatted input for FlyBase scaffolds
bioawk -c fastx '{ print length($seq) }' "$OUTDIR/flybase_scaffolds.fa" \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nFB_Scaff\t0" } { print "FB_Scaff\t" $1 >
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
```

#### Calculate BUSCO scores (compleasm)

I couldn't get BUSCO to work, so I tried compleasm instead. 

First, I made a compleasm env and downloaded it using conda.

```{bash, eval=FALSE}
mamba create -n compleasm_env -c bioconda -c conda-forge compleasm
conda activate compleasm_env
```

Then I made a directory for the compleasm outputs: /code/scripts/out/compleasm.

Finally, I ran compleasm in the following way:

```{bash, eval=FALSE}
compleasm run \
  -a ./out/hifi_fly_assembly.bp.p_ctg.fa \
  -o ./out/compleasm/my_assembly \
  -L diptera \
  -t 16

compleasm run \
  -a ./out/flybase_scaffolds.fa \
  -o ./out/compleasm/flybase_scaffolds \
  -L diptera \
  -t 16
```
Which then stores results inside the summary.txt file inside each of the following paths:
compleasm/my_assembly/
compleasm/flybase_scaffolds/

The BUSCO scores were the following for the scaffold assembly:
```
## lineage: diptera_odb12
S:99.39%, 5035
D:0.47%, 24
F:0.14%, 7
I:0.00%, 0
M:0.00%, 0
N:5066
```

And for my assembly:
```
## lineage: diptera_odb12
S:99.39%, 5035
D:0.39%, 20
F:0.14%, 7
I:0.00%, 0
M:0.08%, 4
N:5066
```
