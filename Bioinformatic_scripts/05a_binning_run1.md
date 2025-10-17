# Bin genomes using aviary workflow

#### First symlink reads from samples to use for differential coverage binning to new folder
##### Samples are grouped based on conditioning tank and treatment, which showed most consistency in sample replicates


```bash
## Tank MIS5, 2Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

for FILE in $(ls $READS_DIR/{SE2451,SE2452,SE2453,SE2454,SE2455,SE2457,SE2473,SE2474,SE2475,SE2476,SE2477,SE2496,SE2497,SE2498,SE2499,SE2500,SE2536,SE2537,SE2538,SE2540,SE2541,SE2542}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/MIS5_2ML_to_bin
done
```


```bash
## Tank MIS7, 2Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

mkdir /srv/projects/microbial_inducers/data/MIS7_2ML_to_bin

for FILE in $(ls $READS_DIR/{SE2458,SE2459,SE2460,SE2461,SE2463,SE2478,SE2479,SE2480,SE2481,SE2482,SE2484,SE2502,SE2503,SE2505,SE2506,SE2507,SE2543,SE2544,SE2545,SE2546,SE2548}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/MIS7_2ML_to_bin
done
```


```bash
## Tank MIS8, 2Month light conditioned

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

mkdir /srv/projects/microbial_inducers/data/MIS8_2ML_to_bin

for FILE in $(ls $READS_DIR/{SE2464,SE2465,SE2468,SE2486,SE2487,SE2488,SE2489,SE2490,SE2491,SE2508,SE2509,SE2510,SE2511,SE2512,SE2513,SE2514,SE2549,SE2550,SE2551,SE2552,SE2553,SE2554}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/MIS8_2ML_to_bin
done
```

## MIS5C assemblies

Running in bactches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS5 (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS5_2ML_to_bin"
# Location of scaffolds/assemblies. Should end in '.fasta'
SCAFFOLD_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
# Location to output results from aviary
AVIARY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning"

mkdir -p $AVIARY_DIR

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

for FILE in $(ls $SCAFFOLD_DIR/{SE2451,SE2454,SE2457,SE2473,SE2477,SE2496}_scaffolds.fasta); do
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

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS5_2ML_to_bin"
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
                        
for FILE in $(ls $SCAFFOLD_DIR/{SE2498,SE2500,SE2536,SE2537,SE2540}_scaffolds.fasta); do
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

Running in bactches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS7 (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS7_2ML_to_bin"
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

for FILE in $(ls $SCAFFOLD_DIR/{SE2458,SE2460,SE2461,SE2478,SE2479,SE2484}_scaffolds.fasta); do
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
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS7_2ML_to_bin"
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

for FILE in $(ls $SCAFFOLD_DIR/{SE2503,SE2506,SE2507,SE2543,SE2544,SE2548}_scaffolds.fasta); do
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

Running in batches as assemblies finish


```bash
#!/bin/bash

# Bin genomes using aviary pipeline

# binning all MIS8 (batch 1)

# load software
module load miniconda3/1.1
conda activate aviary_0.5.0

# Location of short read files
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS8_2ML_to_bin"
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

for FILE in $(ls $SCAFFOLD_DIR/{SE2464,SE2468,SE2485,SE2488,SE2490}_scaffolds.fasta); do
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
SHORT_READ_DIR="/srv/projects/microbial_inducers/data/MIS8_2ML_to_bin"
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

for FILE in $(ls $SCAFFOLD_DIR/{SE2508,SE2510,SE2513,SE2550,SE2553}_scaffolds.fasta); do
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
