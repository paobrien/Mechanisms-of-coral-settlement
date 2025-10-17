# Make protein tree of cycloprodigiosin genes

### Get protein sequences for AGMO annotated genes


```bash
# first get list of genes
grep "K15537" annotations.tsv > K15537_hits.tsv # note: all genes are also annotated as PF04116.16 (the protein family for PRUB680)

# clean header in annotations for exact pattern matching (otherwise matches to incorrect genes)
awk '/^>/ {split($0, a, " "); print a[1]; next} {print}' genes.faa > cleaned_genes.faa

# use list to extract protein seqs (note: grep -Fx for exact line matching; grep -v "^--$" to remove separator line from grep -A)
cut -f1 K15537_hits.tsv | sed 's/^/>/' | grep -Fxf - -A1 cleanded_genes.faa | grep -v "^--$" > K15537_genes.faa
```

### Get PRUB680 protein sequence


```bash
# downloaded from https://www.uniprot.org/uniparc/UPI00026CC750/entry
"/scratch/project/micro_inducers/data/ncbi/cyclo/U1KYF0.fasta"
```

### Get UniProt RP75 alignment for PF04116 

protein family for PRUB680 - AGMO (K15537) also has this annotation


```bash
# alignment downdloaded from https://www.ebi.ac.uk/interpro/entry/pfam/PF04116/ - alignment tab
# includes both full (~20,000 sequences) and seed (curated ~100 sequences) alignments
```

Alignments are in stockholm format, need to convert to fasta


```bash
# convert alignment with HMMER

# reformat
esl-reformat afa PF04116.alignment.full.sto > PF04116.alignment.full.fasta
esl-reformat afa PF04116.alignment.seed.sto > PF04116.alignment.seed.fasta
```

### Add sequences to alignment


```bash
# combine PRUB680 with AGMO seqs
cat K15537_genes.faa U1KYF0.fasta > K15537_U1KYF0.faa
```


```bash
# bunya interactive

# load module
module load mafft/7.490-gcc-10.3.0-with-extensions

# add seqs to alignment
mafft --add K15537_U1KYF0.faa --keeplength --reorder PF04116.alignment.full.fasta > cyclo_combined_full_alignment.fasta
mafft --add K15537_U1KYF0.faa --keeplength --reorder PF04116.alignment.seed.fasta > cyclo_combined_seed_alignment.fasta
```

### Infer phylogenetic tree


```bash
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=36G
#SBATCH --job-name=iqtree
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_ace
#SBATCH -o /home/uqpobri2/scripts/out/iqtree.out
#SBATCH -e /home/uqpobri2/scripts/error/iqtree.error

##################################################

# infer phylogeny using iqtree

# load sofware
module load iqtree/2.2.2.3

# set dir
inDir="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/16_cyclo_tree/alignment"
outDir="/scratch/project/micro_inducers/analysis/ECT03-1_biofilm_MG_2022/16_cyclo_tree/iqtree"

# find best model and infer tree
# --reorder alignment
iqtree -s $inDir/cyclo_combined_seed_alignment.fasta -m TEST -nt 60 -B 1000 --prefix $outDir/cyclo_combined_seed_iqtree
```

### Get tree metadata

Get metadata for protein sequences downloaded from Interpro


```bash
# first get pfam ids from stockholm alignment
grep -v '^#' PF04116.alignment.seed.sto | awk '{print $1}' > pfam_ids.txt
```

Now batch search on the UniProt webpage (https://www.uniprot.org/):

1. Click on 'list' in search bar

2. Upload pfam ids

3. When search is finished click 'Download', change format to TSV.

4. Under UniProt data click the down arrow for 'names and taxonomy'

5. Ensure 'taxonomic lineage' is selected

6. Finally click 'download'




```bash

```
