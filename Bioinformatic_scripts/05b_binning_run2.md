# Bin genomes using aviary workflow

#### First symlink reads from samples to use for differential coverage binning to new folder
##### Samples are grouped based on conditioning tank and treatment, which showed most consistency in sample replicates


```bash
## Tank MIS5, 1Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"

for FILE in $(ls $READS_DIR/{SE6768,SE6769,SE6770,SE6771,SE6772,SE6787,SE6788,SE6789,SE6790,SE6791,SE6792,SE6797,SE6798,SE6799,SE6800,SE6801,SE6802,SE6803,SE6832,SE6833,SE6834}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/MIS5_1ML_to_bin
done
```


```bash
## Tank MIS7, 1Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"

for FILE in $(ls $READS_DIR/{SE6773,SE6774,SE6775,SE6776,SE6777,SE6793,SE6794,SE6795,SE6796,SE6804,SE6805,SE6806,SE6807,SE6808,SE6809,SE6810,SE6820,SE6821,SE6822,SE6835,SE6836,SE6837,SE6838,SE6839,SE6840}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/MIS7_1ML_to_bin
done
```


```bash
## Tank MIS8, 1Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"

for FILE in $(ls $READS_DIR/{SE6778,SE6779,SE6780,SE6781,SE6782,SE6783,SE6784,SE6811,SE6812,SE6813,SE6814,SE6815,SE6816,SE6817,SE6823,SE6824,SE6825,SE6826,SE6827,SE6828,SE6829,SE6841,SE6842,SE6843,SE6844,SE6845,SE6846,SE6847}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/MIS8_1ML_to_bin
done
```


```bash
## Dark 2Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"

for FILE in $(ls $READS_DIR/{SE6848,SE6849,SE6850,SE6851,SE6852,SE6853,SE6854,SE6855,SE6856,SE6857,SE6858,SE6859,SE6860,SE6861,SE6862,SE6863,SE6864,SE6865,SE6866,SE6867}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/dark_to_bin
done
```

## MIS5C assemblies

Running in batches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS5 (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS5_1ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6768,SE6771,SE6788,SE6792,SE6798}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS5 (batch 2)

# load programs
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS5_1ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""
                                
for FILE in $(ls $SCAFFOLD_DIR/{SE6801,SE6803,SE6833,SE6834}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```

## MIS7C assemblies

Running in batches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS7 (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS7_1ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6773,SE6776,SE6793,SE6805,SE6808,SE6810}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS7 (batch 2)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS7_1ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6820,SE6822,SE6836,SE6838,SE6840}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```

## MIS8C assemblies

Running in bactches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS8 (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS8_1ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6778,SE6783,SE6784,SE6811,SE6813,SE6816}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS8 (batch 2)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS8_1ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6825,SE6828,SE6829,SE6842,SE6844,SE6847}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```

## Dark assemblies

Running in batches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all dark (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/dark_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6850,SE6851,SE6853,SE6856,SE6858}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all dark (batch 2)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/dark_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

#mkdir -p $AVIARY_DIR

# Max number of threads to use
THREADS=30
PPTHREADS=30 
# Max memory to allow for bin recovery
MEMORY_RECOVER=300

# Location to store conda environments (don't change)
export CONDA_ENV_PATH=~/.conda/envs/
# Location of GTDBTK database (don't change)
export GTDBTK_DATA_PATH=/srv/db/gtdbtk/official/release207
# Location of EnrichM database (ignore/don't change)
export ENRICHM_DB=""
# Location of EggNog database (ignore/don't change)
export EGGNOG_DATA_DIR=""

for FILE in $(ls $SCAFFOLD_DIR/{SE6860,SE6863,SE6864,SE6867}_scaffolds.fasta); do
#get variables
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "procesessing scaffold file $FILE"
    echo "sample name is $SAMPLE"
# run aviary
    aviary recover \
    --assembly $FILE \
    --paired_reads_1 $SHORT_READ_DIR/*R1.fastq.gz \
    --paired_reads_2 $SHORT_READ_DIR/*R2.fastq.gz \
    --max_threads $THREADS \
    --pplacer_threads $THREADS \
    --max_memory $MEMORY_RECOVER \
    --pplacer_threads $THREADS \
    --n_cores $THREADS \
    --output $AVIARY_DIR/$SAMPLE \
    --workflow recover_mags \
    --semibin-model ocean \
    --conda-frontend "mamba" 2>&1 | tee $AVIARY_DIR/logs/${SAMPLE}_aviary.log
done

```
