# Get bin stats, coverage and taxonomy

### Check percentage of reads that map to bins


```bash
#!/bin/bash

# Map reads to genome bins using coverm

# load software
module load miniconda3
conda activate coverm_0.6.1

READS_DIR="/srv/projects/microbial_inducers/data/all_reads_trimmed"
GENOMES_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95"
OUTPUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/05_coverm/bins_dereplicated_95"

# get genome coverage
coverm genome \
-c $READS_DIR/*fastq.gz \
-d $GENOMES_DIR \
-x fna \
--min-read-aligned-percent 0.75 \
--min-read-percent-identity 0.95 \
-m relative_abundance \
-t 30 \
--output-file $OUTPUT_DIR/dereplicated_genomes_coverage.out \
--output-format dense

```

### Check percentage of scaffolds successfully binned


```bash
#!/bin/bash

# Use singleM to compare single copy markers between assembly and bins

# load software
module load miniconda3/1.1
conda activate singlem_dev_140922

# Using trimmed reads (still sufficient length for singleM)
READS_DIR="/srv/projects/microbial_inducers/data/all_reads_trimmed"
ASSEMBLY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
BINS_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95"
SINGLEM_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/02_singlem/all_samples"

#mkdir -p $SINGLEM_DIR

THREADS=20
METAPACKAGE="/srv/db/singlem_packages/S3.metapackage_20220906.smpkg"

# get single copy markers from bins (better to do this after derep)
singlem pipe \
--metapackage $METAPACKAGE \
--genome-fasta-files $BINS_DIR/*.fna \
--otu-table $SINGLEM_DIR/bins.otu_table.tsv \
--threads $THREADS

singlem appraise \
--metagenome_otu_tables $SINGLEM_DIR/metagenome.otu_table.tsv \
--assembly_otu_tables $SINGLEM_DIR/assembly.otu_table.tsv \
--genome_otu_tables $SINGLEM_DIR/bins.otu_table.tsv \
#--plot $SINGLEM_DIR/singlem_appraise.svg \
&> $SINGLEM_DIR/bins_appraisal.log

```

### Check quality of bins with checkM


```bash
#!/bin/bash

# Get bin quality with checkm

# Load modules
module load miniconda3
conda activate checkm-genome-1.1.3

GENOME_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95"
CHECKM_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/06_checkm/dereplicated_bins_95"

# get checkm quality
checkm lineage_wf \
$GENOME_DIR $CHECKM_DIR -x fna -t 20

# Create checkm file
# short = 1, long = 2

checkm qa \
$CHECKM_DIR/lineage.ms $CHECKM_DIR -o 2 -f $CHECKM_DIR/derep_bins.tsv --tab_table -t 1
```

### Format data files for R input

#### Checkm file


```bash
# replace # with 'No', replace space with '_', remove ()
sed -i -e '1 s/#/No/g;1 s/ /_/g;1 s/[)(]//g' derep_bins.tsv
```

#### Coverm file


```bash
# clean column strings
sed -i -e 's/_trimmed_R1.fastq.gz Relative Abundance (%)//g' dereplicated_genomes_coverage.out
```

### Classify bins using GTDB-tk


```bash
#!/bin/bash

## Get MAG taxonomy with GTDB

# load software
module load miniconda3
conda activate gtdbtk-2.1.0

GENOME_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/07_gtdb/derep_bins_95"

# run gtdb
gtdbtk classify_wf \
--genome_dir $GENOME_DIR \
--out_dir $OUT_DIR \
-x fna \
--prefix derep_bins \
--cpus 20 \
--full_tree

```
