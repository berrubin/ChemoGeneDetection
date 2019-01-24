# README

Script for identifying chemosensory genes in genome sequence.

## Getting Started

### Dependencies

The script runs in python 2.7. You will need biopython installed. You will also need an installation of BLAST and exonerate.

### Using BLAST to find possible hits.

### Running the pipeline

The script ```chemo_finder.py``` launches the pipeline. 

```-p``` specifies the number of cores to use. However, this only has any effect if multiple genomes are being examined. Just a single thread is used for each genome so if two genomes are being examined and two threads are available, then it will run ~twice as fast. Since I developed it to run simultaneously over ~20 genomes this worked perfectly for me but obviously this won't be the case for everyone. Splitting genomes into separate chromosomes and treating them as separate genomes would be an easy way to get around this limitation in lieu of me actually writing a more efficient way to thread it.

```-c``` specifies the full path to the config file which just specifies where the genomes are. There are two tab delimited columns. The first is a four character taxon code and the second is the full path to the genome fasta file.

```-r``` path to the fasta file containing reference protein sequences from the family of genes that you are trying to annotate. 

```-f``` abbreviation for the family of genes being annotated. This isn't actually particularly important as it is just used for labelling purposes. It does need to be consistent however. 

```-i``` specifies the iteration of the pipeline. Again, just needed for labeling purposes.

```-l``` is the minimum length of the protein sequence to be kept as a good sequence. Default is 350. This is generally appropriate for ORs and GRs.

```-s``` alternatively, rather than setting a fixed minimum length, if you're family of genes of interest has a greater diversity of lengths, you can specify a minimum length by the proportion of the reference sequence that it covers. So if you specify ```0.8``` then you would only keep sequences that are at least 80% the length of the reference sequence being used to identify the new sequence. If ```-s``` is specified then ```-l``` is ignored.



```
python chemo_finder.py -p 20 -c genomes.config -f or -i 0 -l 350 -b or_results \
-t tblastn_hits_dir -r all_or_iter_0.fa 
```

### Iterating the pipeline
