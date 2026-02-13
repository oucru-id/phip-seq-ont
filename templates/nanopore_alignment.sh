#!/bin/bash

: '
This aligns reads to the index after trimming the 
read to be the same length as tiles in the library.
'

set -euo pipefail

STREAM_FILE_CMD=!{params.fastq_stream_func}
FASTQ=!{respective_replicate_path}
INDEX=!{index}/peptide
ALIGN_OUT_FN=!{sample_id}.sam
CPUS=!{task.cpus}
OP_ARGS="!{params.bowtie_optional_args}"

$STREAM_FILE_CMD $FASTQ | bowtie2 \
  --threads $CPUS \
  $OP_ARGS \
  -x $INDEX \
  -U - \
  -S $ALIGN_OUT_FN