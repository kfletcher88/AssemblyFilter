A shell wrapper to blast an assembly against a custom database and then positively filter for pre-selected reference sequences.

The wrapper requires the following tools:

1. [NCBI blast](https://www.ncbi.nlm.nih.gov/pubmed/2231712 "Altschul et al. 1990")
2. [Samtools v1.3+](https://www.ncbi.nlm.nih.gov/pubmed/19505943 "Li et al. 2009") 

Input files required are:

1. Query assembly
2. Reference database
3. Pre-selected headers to filter.

For usage type:

```
./AssemblyFilter.sh -h
```

The basic commands called by this shell script, when using NCBI nt, are:

```
blastn -query [Input.fasta] -db [BLASTdb/nt] -max_target_seqs 1 -outfmt 6 -word_size 24 -num_threads 10 -out [Prefix].blastn
awk -v FS='|' '{print $4}' [Prefix].blastn | sed 's/\..*//' | sort -u | comm -12 - <(sort [RefSeqs] | sed 's/\..*//') | join -2 3 - <(sed 's/|/ /3' [Prefix].blastn | sed 's/\..*//' | sort -k3,3) | awk '{print $2}' | sort -u > [Prefix].filt.h
xargs samtools faidx Input.fasta] < [Prefix].filt.h > [Prefix].filt.fasta
```

Where the text contained in square brackets are variables provided to the shell script. Specifically:

`-i` Input.fasta.

`-d` BLAST formatted database, typically nt though a customDB may be supplied if running `-t 2`.

`-l` A list of reference sequences in the database to be positively filtered for.

`-p` Prefix for output files.

Optional arguments include:

`-t` Run type 1 will strip away the gi information from the BLASTn results. Run type 2 will not attempt to do this and is therefore suited for custom DB builds.

`-T` Number of threads for the process to use. The default is set to 10.

`-W` Word size to seed alignments.

The shell script also contains safe-guards to ensure files generated here are not over-written and will clean up intermediate files which are of no use to the user down-stream.
