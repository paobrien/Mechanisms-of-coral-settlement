# Pre-process metagenome data

### Symlink data to projects folder


```bash
# run 1
for FILE in $(ls /srv/data/by_run/NVS031/paulobrien/*/*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/20220921_metagenome_data
done
```


```bash
# run 2
for FILE in $(ls /srv/data/by_run/NVS032/paulobrien/*/*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/20221105_metageome_data_2
done
```


```bash
# run 2 - repeats (linking to same folder)
for FILE in $(ls /srv/data/by_run/NVS046/paulobrien/*/*.fastq.gz); do
    echo $FILE
    ln -s $FILE /srv/projects/microbial_inducers/data/20221105_metageome_data_2
done

```

### Run fastqc to check quality of reads


```bash
#!/bin/bash

# sequencing run 1

# load software
module load miniconda3
conda activate fastqc_0.11.9

# run fastqc for each read file
for FILE in $(ls /srv/data/by_run/NVS031/paulobrien/*/*.fastq.gz); do
    DIR=$(dirname $FILE)
    SAMPLE=$(basename $FILE | cut -d "." -f 1)
# check variables
    echo "filepath is: $DIR"
    echo "fastqc sample: $SAMPLE"
# run fastqc
    fastqc -o /srv/projects/microbial_inducers/analysis/biofilm_metagenome/01_qc_reads/01_fastqc -t 20 --nogroup $DIR/${SAMPLE}.fastq.gz
done

```


```bash
#!/bin/bash

# sequencing run 2

# load software
module load miniconda3
conda activate fastqc_0.11.9

# run fastqc for each read file
for FILE in $(ls /srv/data/by_run/NVS032/paulobrien/*/*.fastq.gz); do
    DIR=$(dirname $FILE)
    SAMPLE=$(basename $FILE | cut -d "." -f 1)
# check variables
    echo "filepath is: $DIR"
    echo "fastqc sample: $SAMPLE"
# run fastqc
  fastqc -o /srv/projects/microbial_inducers/analysis/biofilm_metagenome/01_qc_reads/01_fastqc/run2 -t 20 --nogroup $DIR/${SAMPLE}.fastq.gz
done

```


```bash
#!/bin/bash

# run 2 - repeats

# load software
module load miniconda3
conda activate fastqc_0.11.9

# run fastqc for each read file
for FILE in $(ls /srv/data/by_run/NVS046/paulobrien/*/*.fastq.gz); do
    DIR=$(dirname $FILE)
    SAMPLE=$(basename $FILE | cut -d "." -f 1)
# check variables
    echo "filepath is: $DIR"
    echo "fastqc sample: $SAMPLE"
# run fastqc
  fastqc -o /srv/projects/microbial_inducers/analysis/biofilm_metagenome/01_qc_reads/01_fastqc/run2_repeats -t 20 --nogroup $DIR/${SAMPLE}.fastq.gz
done

```

### Concatenate repeats - some samples sequenced twice due to low data output


```bash
#!/bin/bash

# concatenate sequencing repeats

path="/srv/projects/microbial_inducers/data/20221105_metageome_data_2"
ids="SE6778 SE6779 SE6781 SE6784 SE6789 SE6794 SE6796 SE6810 SE6812 SE6821 SE6824 SE6827 SE6829 SE6835 SE6847 SE6860"

for i in $ids; do
    echo $i
    echo $path/$i
# concatenate file
    zcat $path/${i}_S*_R1_001.fastq.gz > $path/${i}_C000_R1_001.fastq
    zcat $path/${i}_S*_R2_001.fastq.gz > $path/${i}_C000_R2_001.fastq
    gzip $path/${i}_C000_R1_001.fastq
    gzip $path/${i}_C000_R2_001.fastq
done
```

### Quality trim with fastp

note: fastp by default will remove polyG tails from Novaseq artifacts, as well as adapters and quality below Q15


```bash
#!/bin/bash

# run 1

# load software
module load miniconda3
conda activate fastp_0.23.2

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/20220921_metagenome_data"
TRIMMED_DIR="/srv/projects/microbial_inducers/data/20221010_metagenome_data_trimmed"
RESULTS_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/01_qc_reads/03_fastp_trimming"

mkdir $TRIMMED_DIR

# run fastp
for FILE in $(ls $READS_DIR/*R1*.fastq.gz); do
# set variables 
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "trimming sample: $SAMPLE"
# quality trim
    fastp \
    --in1 $READS_DIR/${SAMPLE}_*_R1_*.fastq.gz \
    --in2 $READS_DIR/${SAMPLE}_*_R2_*.fastq.gz \
    --out1 $TRIMMED_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    --out2 $TRIMMED_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -f 15 \
    -F 15 \
    -w 16 \
    -h $RESULTS_DIR/${SAMPLE}.html &> $RESULTS_DIR/${SAMPLE}.log
done

```


```bash
#!/bin/bash

# run 2

# load software
module load miniconda3
conda activate fastp_0.23.2

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/20221105_metageome_data_2"
TRIMMED_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
RESULTS_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/01_qc_reads/03_fastp_trimming/run2"

mkdir $TRIMMED_DIR

# run fastp
for FILE in $(ls $READS_DIR/*R1*.fastq.gz); do
# Set variables. 
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "trimming sample: $SAMPLE"
# quality trim
    fastp \
    --in1 $READS_DIR/${SAMPLE}_*_R1_*.fastq.gz \
    --in2 $READS_DIR/${SAMPLE}_*_R2_*.fastq.gz \
    --out1 $TRIMMED_DIR/${SAMPLE}_trimmed_R1.fastq.gz \
    --out2 $TRIMMED_DIR/${SAMPLE}_trimmed_R2.fastq.gz \
    -f 15 \
    -F 15 \
    -w 16 \
    -h $RESULTS_DIR/${SAMPLE}.html &> $RESULTS_DIR/${SAMPLE}.log
done

```


```bash
#!/bin/bash

# run 2 - repeats merged

# load software
module load miniconda3
conda activate fastp_0.23.2

# set directories
READS_DIR="/srv/projects/microbial_inducers/data/20221105_metageome_data_2"
TRIMMED_DIR="/srv/projects/microbial_inducers/data/20221206_metagenome_data_2_trimmed"
RESULTS_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/01_qc_reads/03_fastp_trimming/run2"

#mkdir $TRIMMED_DIR

# run fastp
for FILE in $(ls $READS_DIR/*_C000_R1*.fastq.gz); do
# Set variables. 
    SAMPLE=$(basename $FILE | cut -d "_" -f 1)
    echo "trimming sample: $SAMPLE"
# quality trim
    fastp \
    --in1 $READS_DIR/${SAMPLE}_R1_*.fastq.gz \
    --in2 $READS_DIR/${SAMPLE}_R2_*.fastq.gz \
    --out1 $TRIMMED_DIR/${SAMPLE}_trimmedc000_R1.fastq.gz \
    --out2 $TRIMMED_DIR/${SAMPLE}_trimmedc000_R2.fastq.gz \
    -f 15 \
    -F 15 \
    -w 16 \
    -h $RESULTS_DIR/${SAMPLE}.html &> $RESULTS_DIR/${SAMPLE}.log
done

```
