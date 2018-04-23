#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to final directory
novo_final=${project_dir}/novo_final

samtools merge ${novo_final}/MGD3_SO_CAGATC_novo_merge_novo_final.bam \
${novo_final}/MGD2_SO_CAGATC_novo_merge_novo_final.bam \
${novo_final}/MGD_SO_CAGATC_novo_merge_novo_final.bam

#move unmerged away
mkdir ${novo_final}/Anc_unmerged
mv ${novo_final}/ MGD2_SO_CAGATC_novo_merge_novo_final.bam ${novo_final}/Anc_unmerged
mv ${novo_final}/ MGD_SO_CAGATC_novo_merge_novo_final.bam ${novo_final}/Anc_unmerged
```
