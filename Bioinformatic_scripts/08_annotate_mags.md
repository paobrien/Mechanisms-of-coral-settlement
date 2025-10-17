# Annotate MAGs with DRAM 


```bash
#!/bin bash

# load program
module load miniconda3
conda activate dram_1.3.3

# set variables
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/08_annotate/dram"
GTDB_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/07_gtdb/derep_bins_95/classify"
CHECKM_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/06_checkm/dereplicated_bins_95"

#note: genome_dir needs to be written directly due '' syntax

# run annotation pipeline
DRAM.py annotate \
-i '/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95/*.fna' \ 
-o $OUT_DIR/dram_annotate \
--gtdb_taxonomy $GTDB_DIR/derep_bins.bac120.summary.tsv \
--gtdb_taxonomy $GTDB_DIR/derep_bins.ar53.summary.tsv \
--checkm_quality $CHECKM_DIR/derep_bins.tsv \
--threads 20 \
--verbose

# run distill to summarise annotation
DRAM.py distill \
-i $OUT_DIR/dram_annotate/annotations.tsv \
-o $OUT_DIR/dram_distill \
--rrna_path $OUT_DIR/dram_annotate/rrnas.tsv \
--trna_path $OUT_DIR/dram_annotate/trnas.tsv \
--groupby_column fasta

```
