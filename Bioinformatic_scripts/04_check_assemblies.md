# Check assembly quality

### Map reads to assemblies to check proportion assembled

note: results are for each sample back to each individual assembly


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate coverm_0.6.1

# set directories
SCAFFOLDS_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
READS_DIR="/srv/projects/microbial_inducers/data/all_reads_trimmed"
OUTPUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/05_coverm/assemblies"

# get % reads mapped to assemblies
for FILE in $(ls $SCAFFOLDS_DIR/*fasta); do
    echo "file is $FILE"
    SAMPLE=$(basename $FILE | cut -d '_' -f 1)
    echo "sample is $SAMPLE"
# run coverm
    coverm genome \
    -c $READS_DIR/$SAMPLE*.fastq.gz \
    -f $FILE \
    -x fasta \
    -t 30 \
    -p minimap2-sr \
    -m relative_abundance \
    --min-read-aligned-percent 0.75 \
    --min-read-percent-identity 0.95 \
    --output-file $OUTPUT_DIR/${SAMPLE}_coverm_contig.out \
    --output-format dense
done
```

### Check proportion of single copy marker genes in assembly


```bash
#!/bin/bash

# load software
module load miniconda3/1.1
conda activate singlem

# using trimmed reads (still sufficient length for singleM)
READS_DIR="/srv/projects/microbial_inducers/data/all_reads_trimmed"
ASSEMBLY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
SINGLEM_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/02_singlem/all_samples"

mkdir -p $SINGLEM_DIR

THREADS=20
METAPACKAGE="/srv/db/singlem_packages/S3.metapackage_20220906.smpkg"

# get single copy markers from assemblies
singlem pipe \
--metapackage $METAPACKAGE \
--genome-fasta-files $ASSEMBLY_DIR/*.fasta \
--otu-table $SINGLEM_DIR/assembly.otu_table.tsv \
--threads $THREADS

# compare against reads
singlem appraise \
--metagenome-otu-tables $SINGLEM_DIR/metagenome.otu_table.tsv \
--assembly-otu-tables $SINGLEM_DIR/assembly.otu_table.tsv \
#--genome-otu-tables $SINGLEM_DIR/bins.otu_table.tsv \
#--plot $SINGLEM_DIR/singlem_appraise.svg \
#--output-assembled-otu-table $SINGLEM_DIR/assembled.otu_table.tsv \
&> $SINGLEM_DIR/assembly_appraisal.log

```
