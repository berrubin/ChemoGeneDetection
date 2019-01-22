# README

Script for identifying chemosensory genes in genome sequence.

## Getting Started

### Dependencies

The script runs in python 2.7. You will need biopython installed. You will also need an installation of BLAST and exonerate.

### Using BLAST to find possible hits.

### Running the pipeline

The rest of the python scripts synthesize these pairwise genome alignments into multiple genome alignments across taxa. ```maf_parser.py``` launches the pipeline. We used the command:

```
python maf_parser.py -i alignments -s AMEL -t AFLO,EMEX,BIMP,BTER,MQUA,HLAB,CCAL,MROT,LALB,DNOV -f gffs \
-w 500 -z 250 -o AMEL_base -m 8 -c bees.tree -p 14 -k AMEL_all_scafs.txt
```

### Iterating the pipeline
