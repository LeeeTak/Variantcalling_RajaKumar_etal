#!/bin/bash

threadn=$1

fastqc *trimmed.fq.gz --threads $threadn --outdir fastqc_trimmed
