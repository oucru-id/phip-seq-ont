#!/bin/bash

set -euo pipefail
CPUS=!{task.cpus}

samtools view -u -@ ${CPUS} !{sam_file} | \
    samtools sort -@ ${CPUS} - > !{sample_id}.bam

samtools sort -@ ${CPUS} !{sample_id}.bam -o !{sample_id}.sorted 

mv !{sample_id}.sorted !{sample_id}.bam

samtools index -b !{sample_id}.bam

samtools idxstats !{sample_id}.bam | \
    cut -f 1,3 | \
    sed "/^*/d" > !{sample_id}.counts
