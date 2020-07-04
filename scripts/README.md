## Running the stacks pipeline
The `run_stacks.sh` script is primarily for documentation purposes and contains no exception handling. User arguments: 1) A file listing all input fastq files with one file per line, 2) a trimmomatic adapter fasta file, 3) number of cores to use for stacks, 4) the enzyme recognition site to check at the start of each read, and 4) maximum stacks distance (-n, -M).
```
./run_stacks.sh fastq_list.txt adapters.fa 24 CGT 3
```
