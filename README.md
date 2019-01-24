# README

Script for identifying chemosensory genes in genome sequence. The inspiration and steps are largely drawn from:

Zhou et al. 2012 Phylogenetic and transcriptomic analysis of chemosensory receptors in a pair of divergent ant species reveals sex-specific signatures of odor coding. PLoS Genet. 8:e1002930.

and

Zhou et al. 2015 Chemoreceptor evolution in Hymenoptera and its implications for the evolution of eusociality. GBE. 7(8): 2407–2416.

However, there are many iterations of this type of pipeline and all seem generally effective. I don't expect that this version is substantially better (and hopefully isn't much worse) than those used previously.

## Getting Started

### Dependencies

The script runs in ```python 2.7```. You will need ```biopython``` installed. You will also need an installation of ```BLAST``` and ```exonerate```.

### Using BLAST to find possible hits.
The first step is to use ```TBLASTN``` to find regions that potentially contain a full length gene of interest. There is nothing special about these commands but the output does need to be in a particular format (```-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe"```). So, for example:

```makeblastdb -in AMEL_genome.fasta -out AMEL_genome_db -dbtype nucl```

followed by:

```tblastn -query or_ref_seqs.fa -db AMEL_genome_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe" -out tblastn_hits_dir/AMEL_or_iter_0_1000.txt -evalue 1000```

Where ```AMEL_genome.fasta``` is a fasta file of the genome of interest and ```or_ref_seqs.fa``` is fasta formatted protein sequences from the gene family of interest.

### Running the pipeline

The script ```chemo_finder.py``` launches the pipeline. 

#### Required parameters
```-c``` specifies the full path to the config file which just specifies where the genomes are. There are two tab delimited columns. The first is a four character taxon code and the second is the full path to the genome fasta file.

```-r``` path to the fasta file containing reference protein sequences from the family of genes that you are trying to annotate. 

```-b``` base directory where output will be written.

```-t``` directory that contains results from ```TBLASTN``` for current iteration.

```-f``` abbreviation for the family of genes being annotated. This isn't actually particularly important as it is just used for labelling purposes. It does need to be consistent however. 

```-i``` specifies the iteration of the pipeline. Again, just needed for labeling purposes.

```-x``` path to directory containing exonerate executable. If it is already in your path then there is no need to specify here.

#### Optional parameters
```-p``` specifies the number of cores to use. However, this only has any effect if multiple genomes are being examined. Just a single thread is used for each genome so if two genomes are being examined and two threads are available, then it will run ~twice as fast. Since I developed it to run simultaneously over ~20 genomes this worked perfectly for me but obviously this won't be the case for everyone. Splitting genomes into separate chromosomes and treating them as separate genomes would be an easy way to get around this limitation in lieu of me actually writing a more efficient way to thread it.

```-l``` is the minimum length of the protein sequence to be kept as a good sequence. Default is 350. This is generally appropriate for ORs and GRs.

```-s``` alternatively, rather than setting a fixed minimum length, if you're family of genes of interest has a greater diversity of lengths, you can specify a minimum length by the proportion of the reference sequence that it covers. So if you specify ```0.8``` then you would only keep sequences that are at least 80% the length of the reference sequence being used to identify the new sequence. If ```-s``` is specified then ```-l``` is ignored.

```-k``` the pipeline gets seed locations for identifying genes of interest using ```TBLASTN``` which might miss more divergent matches. This parameter allows you to include flanking sequence adjacent to the ```BLAST``` hits. Then this flanking sequence will be included when the more thorough annotation is done with ```exonerate```. The bigger this is the slower it will run but the less likely you are to miss the ends of genes. Default is 1,000.

```-n``` maximum intron size. Again, the bigger this is, the longer it will take for ```exonerate``` alignments but the less likely you are to miss genes. Default is 12,000.

```-e``` maximum evalue of ```BLAST``` hits to examine. This should be pretty high, especially for the early iterations where you are basing annotations on highly divergent proteins. Default is 100.

```-m``` when choosing between possible annotations of a particular gene, we often want the longer one, regardless of the similarity score between the reference sequence and the sequence of interest. However, sometimes you want to take the more similar sequence instead, particularly when lengths are similar. When two possible annotations differ by this length or less, then the one with the higher score is used. Default is 40.

```-v``` the amount of overlap when merging query sequence based on HSPs. Default is 30. Probably don't want to change this.

Example command:
```
python chemo_finder.py -p 2 -c genomes.config -f or -i 0 -l 350 -b or_results \
-t tblastn_hits_dir -r or_ref_seqs.fa 
```

### Output
A directory will be created with the name specified by ```-b```. Within this directory will be two other directories. The names of these are ```[family]_genes/regional_iter_[iteration]_[evalue]_[flank size]_[max intron size]_[max overlap]_[allowed length differential]```. The ```regional``` directory contains that results for individual genes from ```exonerate```. The ```genes``` directory contains the main results. For each species, there are four files containing the CDS ("fna"), the protein sequence ("faa"), the transcript sequence ("trans.fna"), and the genomic coordinates ("gff") for all of the identified genes. There is also a directory for each species with the reconstructed sequences for each annotated gene which can generally be ignored.

### Iterating the pipeline
Running the pipeline once is not going to recover all of the genes of interest. Searching for signatures of the sequences recovered in the first iteration can yield many more genes. This is especially useful if you are examining many genomes simultaneously so that they can inform each other. Therefore, I would recommend combining all of the proteins identified in each iteration (plus the original reference sequences) and starting the process over again by BLASTing these annotated proteins against the genome. For odorant receptors in halictids, I found that three iterations was generally optimal but this will likely vary.

### Notes
If anybody wants this to work with Genewise, let me know. I have the code for it but it didn't work as well for me as exonerate did so I didn't include it here.