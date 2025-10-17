# Assemble metagenome using MetaSpades 

#### Symlink samples to assemble to new directory

Assembling one replicate from each Treatment x tank x settlement combination first to save processing time


```bash
# Psin samples

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

for FILE in $(ls $READS_DIR/{SE2451,SE2454,SE2457,SE2458,SE2460,SE2461,SE2464,SE2468}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/psin_to_assemble
done
```


```bash
# Dfav

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

for FILE in $(ls $READS_DIR/{SE2473,SE2477,SE2478,SE2479,SE2484,SE2485,SE2488,SE2490}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/dfav_to_assemble
done
```


```bash
# Plob

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

for FILE in $(ls $READS_DIR/{SE2496,SE2498,SE2500,SE2503,SE2506,SE2507,SE2508,SE2510,SE2513}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/plob_to_assemble
done
```


```bash
# Easp

READS_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"

for FILE in $(ls $READS_DIR/{SE2536,SE2537,SE2540,SE2543,SE2544,SE2548,SE2550,SE2553}*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/easp_to_assemble
done
```

### Assemble metagenomes - sequencing run1

#### Psin samples


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/psin_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/psin"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 600
done

```

#### Dfav


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/dfav_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/dfav"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 300
done

```


```bash
#!/bin/bash

# re-running SE2485 - out of memory

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/dfav_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/dfav"

# run metaspades
time metaspades.py -o $OUT_DIR/SE2485 -t 20 -m 600 --restart-from last

conda deactivate

```

#### Plob


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/plob_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/plob"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
# get variables
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    nice -n 10 metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 400
done

```

#### Easp


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/easp_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/easp"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
# get variables
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 550
done

```


```bash
#!/bin/bash

# re-running SE2553 as ran out of memory

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/easp_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/easp"

# run metaspades
time metaspades.py -o $OUT_DIR/SE2553 -t 20 -m 750 --restart-from last

```

### Consolidate scaffolds and compress assembly files

#### Psin


```bash
for DIR in $(ls | grep -v "tar"); do
    echo $DIR
# get file
    ls $DIR/scaffolds.fasta
    ls $DIR/spades.log
# mv file
    cp $DIR/scaffolds.fasta ../all_scaffolds/${DIR}_scaffolds.fasta
    cp $DIR/spades.log ../all_logs/${DIR}_spades.log
done

# compress
for DIR in $(ls | grep -v "tar"); do
    echo "compressing $DIR"
    tar -zcvf ${DIR}.tar.gz $DIR --remove-files $DIR
done
```

#### Dfav


```bash
for DIR in $(ls | grep -v "tar"); do
    echo $DIR
# get file
    ls $DIR/scaffolds.fasta
    ls $DIR/spades.log
# mv file
    cp $DIR/scaffolds.fasta ../all_scaffolds/${DIR}_scaffolds.fasta
    cp $DIR/spades.log ../all_logs/${DIR}_spades.log
# compress
    echo "compressing $DIR"
    tar -zcvf ${DIR}.tar.gz $DIR --remove-files $DIR
done
```

#### Plob


```bash
for DIR in $(ls | grep -v "tar"); do
    echo $DIR
# get file
    ls $DIR/scaffolds.fasta
    ls $DIR/spades.log
# mv file
    cp $DIR/scaffolds.fasta ../all_scaffolds/${DIR}_scaffolds.fasta
    cp $DIR/spades.log ../all_logs/${DIR}_spades.log
# compress
    echo "compressing $DIR"
    tar -zcvf ${DIR}.tar.gz $DIR --remove-files $DIR
done
```

#### Easp


```bash
for DIR in $(ls | grep -v "tar"); do
# check variables
    echo $DIR
# get file
    ls $DIR/scaffolds.fasta
    ls $DIR/spades.log
# mv file
    cp $DIR/scaffolds.fasta ../all_scaffolds/${DIR}_scaffolds.fasta
    cp $DIR/spades.log ../all_logs/${DIR}_spades.log
# compress
    echo "compressing $DIR"
    tar -zcvf ${DIR}.tar.gz $DIR --remove-files $DIR
done
```
