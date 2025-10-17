# Assemble metagenome using MetaSpades 

#### Symlink samples to new directory


```bash
# Psin samples

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
LINK_DIR="/srv/projects/microbial_inducers/data/psin_to_assemble2"

mkdir $LINK_DIR

for FILE in $(ls $READS_DIR/{SE6798,SE6801,SE6803,SE6805,SE6808,SE6810,SE6811,SE6813,SE6816}*.fastq.gz); do
    echo $FILE
    ln -s $FILE $LINK_DIR
done
```


```bash
# Dfav

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
LINK_DIR="/srv/projects/microbial_inducers/data/dfav_to_assemble2"

mkdir $LINK_DIR

for FILE in $(ls $READS_DIR/{SE6768,SE6771,SE6773,SE6776,SE6778,SE6783,SE6784}*.fastq.gz); do
    echo $FILE
    ln -s $FILE $LINK_DIR
done
```


```bash
# Plob

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
LINK_DIR="/srv/projects/microbial_inducers/data/plob_to_assemble2"

mkdir $LINK_DIR

for FILE in $(ls $READS_DIR/{SE6788,SE6792,SE6793,SE6820,SE6822,SE6825,SE6828,SE6829}*.fastq.gz); do
    echo $FILE
    ln -s $FILE $LINK_DIR
done
```


```bash
# Easp

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
LINK_DIR="/srv/projects/microbial_inducers/data/easp_to_assemble2"

mkdir $LINK_DIR

for FILE in $(ls $READS_DIR/{SE6833,SE6834,SE6836,SE6838,SE6840,SE6842,SE6844,SE6847}*.fastq.gz); do
    echo $FILE
    ln -s $FILE $LINK_DIR
done
```


```bash
# Dark biofilms

READS_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
LINK_DIR="/srv/projects/microbial_inducers/data/dark_to_assemble"

mkdir $LINK_DIR

for FILE in $(ls $READS_DIR/{SE6850,SE6851,SE6853,SE6856,SE6858,SE6860,SE6863,SE6864,SE6867}*.fastq.gz); do
    echo $FILE
    ln -s $FILE $LINK_DIR
done
```

### Assemble metagenomes - sequencing run2

#### Psin samples


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/psin_to_assemble2"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/psin/run2"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 500
done

```

#### Dfav


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/dfav_to_assemble2"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/dfav/run2"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 500
done

```

#### Plob


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/plob_to_assemble2"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/plob/run2"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 500
done

```

#### Easp


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/easp_to_assemble2"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/easp/run2"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 500
done

```

#### Dark


```bash
#!/bin/bash

# load software
module load miniconda3
conda activate spades_3.15.3

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/dark_to_assemble"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/dark"

for FILE in $(ls $READS_DIR/*R1.fastq.gz); do
    echo $FILE
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo $SAMPLE
# run metaspades
    time metaspades.py \
    -1 $READS_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    -2 $READS_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -o $OUT_DIR/${SAMPLE} \
    -t 20 -m 500
done

```

### Consolidate scaffolds and compress assembly files


```bash
for DIR in $(ls | grep -v "tar"); do
    echo $DIR
# get file
    ls $DIR/scaffolds.fasta
    ls $DIR/spades.log
# mv file
    cp $DIR/scaffolds.fasta ../../all_scaffolds/${DIR}_scaffolds.fasta
    cp $DIR/spades.log ../../all_logs/${DIR}_spades.log
done

# compress
for DIR in $(ls | grep -v "tar"); do
    echo "compressing $DIR"
    tar -zcvf ${DIR}.tar.gz $DIR --remove-files $DIR
done

```

For darks


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
