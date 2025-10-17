# Taxonomic profiling of reads

##### Use singleM to get a taxonomic profile of raw reads. Note: singleM suggests to use untrimmed reads to keep sequence length. However, using trimmed reads since only ~15bp removed (sequencing noise)

### Symlink reads from both runs to new directory


```bash
# run 1
data_dir="/srv/projects/microbial_inducers/data/"

for FILE in $(ls $data_dir/20221010_metagenome_data_trimmed/*.fastq.gz); do
    echo $FILE
    ln -s $FILE $data_dir/all_reads_trimmed
done

# run2
for FILE in $(ls $data_dir/20221206_metagenome_data_2_trimmed/*.fastq.gz); do
    echo $FILE
    ln -s $FILE $data_dir/all_reads_trimmed
done
```

#### Remove unconcatenated repeat files


```bash
data_dir="/srv/projects/microbial_inducers/data/all_reads_trimmed"

for FILE in $(ls $data_dir/*c000_R1.fastq.gz); do
    echo $FILE
    ID=$(basename $FILE | cut -d '_' -f 1)
    echo $ID
    rm $data_dir/${ID}_trimmed_R1.fastq.gz
    rm $data_dir/${ID}_trimmed_R2.fastq.gz
done
```

### Run singleM


```bash
#!/bin/bash

# load software
module load miniconda3/1.1
conda activate singlem

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/all_reads_trimmed"
SINGLEM_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/02_singlem"

# set variables
THREADS=20
METAPACKAGE="/srv/db/singlem_packages/S3.metapackage_20220906.smpkg"


# generate OTU table using forward and reverse reads for all samples
# note: --archive-otu-table needed as input for condense command
singlem pipe \
--singlem-metapackage $METAPACKAGE \
--forward $READS_DIR/*R1.fastq.gz \
--reverse $READS_DIR/*R2.fastq.gz \
--otu-table $SINGLEM_DIR/all_samples/metagenome.otu_table.tsv \
--archive-otu-table $SINGLEM_DIR/all_samples/archive.otu_table.tsv \
--threads $THREADS \
--assignment-threads $THREADS

# combine OTU tables across multiple markers into a single OTU table.
singlem condense \
--singlem-metapackage $METAPACKAGE \
--input-archive-otu-tables $SINGLEM_DIR/all_samples/archive.otu_table.tsv \
--output-otu-table $SINGLEM_DIR/all_samples/metagenome_condensed.tsv \
--krona $SINGLEM_DIR/metagenome_krona

```
