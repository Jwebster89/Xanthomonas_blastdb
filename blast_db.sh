#!/bin/bash
f=$(basename $1)
blastn -db ./XanthDB -query $1 -outfmt "7 qacc bitscore length pident sacc stitle" -num_threads 8 -max_target_seqs 5 -max_hsps 5  > blast_output/${f}_blast_results.txt
