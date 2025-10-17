# Blast gene sequences to MAGs 

##### Note: make a database of assemblies and MAGs, since sequences may not have been binned


```bash
#!/bin bash

# load software
module load miniconda3
conda activate blast_2.12.0 

MAG_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/04_binnning/dereplicated_bins_95"
ASSEMBLY_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/03_assembly/all_scaffolds"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/09_blast/blast_db"

## make blast database 

# copy MAGS to new directory and append MAG name to contigs
for FILE in $(ls $MAG_DIR/*.fna); do
    # get file name
    NAME=$(basename -s .fna $FILE)
    echo "file name is $NAME"
   # copy file to out dir and append contig names
    cp $FILE $OUT_DIR
    sed -i -e "s/^>/>${NAME}_/" $OUT_DIR/${NAME}.fna
done

# concatenate to one file
cat $OUT_DIR/*.fna > $OUT_DIR/all_mags.fasta

# remove copied mags
for FILE in $(ls $OUT_DIR/*.fna); do
    echo $FILE
    rm $FILE
done

# make MAG database
makeblastdb \
    -in $OUT_DIR/all_mags.fasta \
    -dbtype nucl \
    -out $OUT_DIR/all_mags_db


## create assembly database

# concatenate to one file
cat $ASSEMBLY_DIR/*.fasta > $OUT_DIR/all_scaffolds.fasta

# make scaffold database
makeblastdb \
    -in $OUT_DIR/all_scaffolds.fasta \
    -dbtype nucl \
    -out $OUT_DIR/all_scaffolds_db
```

## Run blast on TBP biosynthesis genes

Using tblastx as per Alker et al. 2020 - use translated proteins rather than nucleotides

Considers all possible codon translations, accounting for the possibility of different reading frames and possible frameshifts.


```bash
#!/bin bash

# run blast

# load software
module load miniconda3
conda activate blast_2.12.0 

DB_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/09_blast/blast_db"
OUT_DIR="/srv/projects/microbial_inducers/analysis/biofilm_metagenome/09_blast/blast_results/tbp_biosynthesis"
DATA_DIR="/srv/projects/microbial_inducers/data/ncbi"

# blast to MAG database
tblastx \
-query $DATA_DIR/pseudoalteromonas_bmp_genes.txt \
-db $DB_DIR/all_mags/all_mags_db \
-out $OUT_DIR/tblastx_results_tbp_to_mags.txt \
-outfmt 6

# blast to scaffold database
tblastx \
-query $DATA_DIR/pseudoalteromonas_bmp_genes.txt \
-db $DB_DIR/all_scaffolds/all_scaffolds_db \
-out $OUT_DIR/tblastx_results_tbp_to_scaffolds.txt \
-outfmt 6

conda deactivate
```

### Find which genomes have a signficant match to TBP genes


```bash
# output significant results, in this case percent identity > 30, aligment >100, e-value < 1e-03, bit score >50
# based on Pearson 2013

# get significant results
awk '$3 > 30 && $4 > 100 && $11 < 1e-03 && $12 > 50 {print}' tblastx_results_tbp_to_mags.txt | \
sort -k 12,12nr > tblastx_significant_results_tbp_to_mags_e03-100.txt

# get unique gene/MAG matches
awk '{print $1, $2}' tblastx_significant_results_tbp_to_mags_e03-100.txt | sort | uniq > tblastx_unique_matches_tbp_to_mags_e03-100.txt

```

## Run blast on cycloprodigiosin gene

##### Running on Bunya server (decommisioning ACE servers)

#### Blast Cyclo nucleotide sequence (downloaded from JBEI, de Rond et al. 2017) 


```bash
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --job-name=blast
#SBATCH --time=4:00:00
#SBATCH --partition=general
#SBATCH --account=a_ace
#SBATCH -o /home/uqpobri2/scripts/out/blast.output
#SBATCH -e /home/uqpobri2/scripts/error/blast.error

##########################

# load software
module load blast

# assign variables
dbDir="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/09_blast/db"
DataDir="/scratch/project/micro_inducers/data/ncbi/cyclo"
OutDir="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/09_blast/blast_results/cyclo"

mkdir -p $OutDir

# check variables
echo ""
echo "#====================#"
echo "#===== JOB INFO =====#"
echo "#====================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "    Input genome path :  ${dbDir}"
echo "    Input sequences   :  ${DataDir}" 
echo "    Output directory  :  ${OutDir}"
echo ""
echo "#====================#"
echo ""

# blast to MAG database
srun tblastx \
    -query $DataDir/PRUB680_cycloprodigiosin.fasta \
    -db $dbDir/all_mags_db \
    -out $OutDir/tblastx_results_cyclo_to_mags.txt \
    -outfmt 6

# blast to scaffold database
srun tblastx \
    -query $DataDir/PRUB680_cycloprodigiosin.fasta \
    -db $dbDir/all_scaffolds_db \
    -out $OutDir/tblastx_results_cyclo_to_scaffolds.txt \
    -outfmt 6

```

#### Find which genomes have a signficant match to cyclo genes


```bash
# output significant results, in this case percent identity > 30, aligment >100, e-value < 1e-06, bit score >50
# based on Pearson 2013

# get significant results
awk '$3 > 30 && $4 > 100 && $11 < 1e-06 && $12 > 50 {print}' tblastx_results_cyclo_to_mags.txt | \
sort -k 12,12nr > tblastx_significant_results_cyclo_to_mags_e06-100.txt

# get unique gene/MAG matches
awk '{print $1, $2}' tblastx_significant_results_cyclo_to_mags_e03-100.txt | sort | uniq > tblastx_unique_matches_cycl_to_mags_e06-100.txt

```

#### Repeat with full biosynthesis pathway (pigA-J, downloaded from KEGG) 


```bash
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --job-name=blast
#SBATCH --time=4:00:00
#SBATCH --partition=general
#SBATCH --account=a_ace
#SBATCH -o /home/uqpobri2/scripts/out/blast.output
#SBATCH -e /home/uqpobri2/scripts/error/blast.error

##########################

# load software
module load blast

# assign variables
dbDir="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/09_blast/db"
DataDir="/scratch/project/micro_inducers/data/ncbi/cyclo"
OutDir="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/09_blast/blast_results/cyclo"

mkdir -p $OutDir

# check variables
echo ""
echo "#====================#"
echo "#===== JOB INFO =====#"
echo "#====================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "    Input genome path :  ${dbDir}"
echo "    Input sequences   :  ${DataDir}" 
echo "    Output directory  :  ${OutDir}"
echo ""
echo "#====================#"
echo ""

# blast to MAG database
srun tblastx \
    -query $DataDir/cyclo_full_short_names.fasta  \
    -db $dbDir/all_mags_db \
    -out $OutDir/tblastx_results_cyclo_full_to_mags.txt \
    -outfmt 6

# blast to scaffold database
srun tblastx \
    -query $DataDir/cyclo_full_short_names.fasta  \
    -db $dbDir/all_scaffolds_db \
    -out $OutDir/tblastx_results_cyclo_full_to_scaffolds.txt \
    -outfmt 6

```

#### Find which genomes have a signficant matches


```bash
# output significant results, in this case percent identity > 30, aligment >100, e-value < 1e-06, bit score >50
# based on Pearson 2013

# get significant results
awk '$3 > 30 && $4 > 100 && $11 < 1e-06 && $12 > 50 {print}' tblastx_results_cyclo_full_to_mags.txt | \
sort -k 12,12nr > tblastx_significant_results_cyclo_full_to_mags_e06-100.txt

```
