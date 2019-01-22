# README

Script for identifying chemosensory genes in genome sequence.

## Getting Started

### Dependencies

The script runs in python 2.7. You will need biopython installed. You will also need an installation of BLAST and exonerate.

### Using BLAST to find possible hits.

### Running the pipeline

The script ```chemo_finder.py``` launches the pipeline. 

```
python chemo_finder.py -p 20 -c genomes.config -f or -i 0 -l 350 -b or_results \
-t tblastn_hits_dir -r all_or_iter_0.fa 
```

### Iterating the pipeline
