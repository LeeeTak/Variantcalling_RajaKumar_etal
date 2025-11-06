#!/usr/bin/bash

files=$(paste -sd "," fqfiles)
Trinity --seqType fq --single $files --max_memory 50G --CPU 8
