#!/usr/bin/bash

funannotate predict -i $1 \
    --species "Plectosphaerella cucumerina" \
    --isolate Plecuc1 \
    --transcript_evidence Trinity/trinity_out_dir/Trinity.fasta \
    -o funannotate
