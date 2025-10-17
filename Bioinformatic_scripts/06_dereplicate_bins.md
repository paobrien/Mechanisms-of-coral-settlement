# Consolidate and dereplicate bins from multiple samples

### Consolodate bins from sample replicates


```bash
#!/bin/bash

# sym link and rename bins to one directory

AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"
LINK_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/all_bins"

mkdir $LINK_DIR

for sample_dir in $(ls $AVIARY_DIR | grep -E -v 'all_bins|_test|logs'); do
    echo $sample_dir
    if [[ -d $AVIARY_DIR/$sample_dir ]]; then
        for bin_file in $(ls $AVIARY_DIR/$sample_dir/bins/final_bins/*.fna); do
            echo $bin_file
            bin_file_basename=$(basename $bin_file .fna | sed "s/.tsv//; s/_bins//")
            echo $bin_file_basename
            sample_out_name=$(echo ${sample_dir}.${bin_file_basename})
            echo $sample_out_name
            ln -s $bin_file $LINK_DIR/${sample_out_name}.fna
        done
    fi
done

```

## Dereplicate bins


```bash
#!/bin/bash

# Dereplicate genomes at 95ANI using coverm

# load software
module load miniconda3
conda activate checkm-genome-1.1.3

GENOME_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/all_bins"
CHECKM_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/06_checkm"
DEREP_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95"

## first need checkm results ##

# get checkm quality
checkm lineage_wf \
$GENOME_DIR $CHECKM_DIR -x fna -t 32

# Create checkm file
# Note: for use in coverm cluster (dereplication), need short checkm output
# short = 1, long = 2

checkm qa \
$CHECKM_DIR/lineage.ms $CHECKM_DIR -o 1 -f $CHECKM_DIR/all_bins.tsv --tab_table -t 1

## dereplicate bins ##

# load software
conda activate coverm_0.6.1

# dereplicate with coverm
coverm cluster \
--genome-fasta-directory $GENOME_DIR \
-x fna \
--ani 95 \
--checkm-tab-table $CHECKM_DIR/all_bins.tsv \
--output-representative-fasta-directory $DEREP_DIR \
--precluster-method finch \
--min-completeness 50 \
--max-contamination 10 \
-t 32

```
