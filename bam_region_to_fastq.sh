#!/usr/bin/env bash

module load htslib/1.9
module load samtools/1.9



#
##  Globals
#
n=1000000
out=fastq_region


mkdir -p $out

region=Chr09:35542701-35551878


#
## Extract Guppy Fastq reads
#
# guppy_fastq=/workspace/hrpelg/Red_Flesh_ON/Guppy_basecalling/02.poreChop/After_PoreCHOP_RedFlesh_ON_GUPPY_cas.fastq.gz
fastq=fastq/After_PoreCHOP_RedFlesh_ON_GUPPY_cas.fastq.gz
  bam=/workspace/hrpelg/Red_Flesh_ON/Guppy_basecalling/04.canu/minimap/correct_reads_to_GD/Guppy_canu_corrected_vs_GDv1.1.bam

quids=$(samtools view $bam $region | perl -lane 'print $F[0]' | xargs)
samtools fqidx -n $n $fastq $quids > >(bgzip > $out/$(basename $fastq))




#
## Extract Albacore Fastq reads
#
fastq=fastq/All_DS_RedFlesh_ON_run1_cas_after_porechop_dis.fastq.gz
  bam=/workspace/hrpelg/Red_Flesh_ON/03.minimap2/Redflesh.bam

quids=$(samtools view $bam $region | perl -lane 'print $F[0]' | xargs)
samtools fqidx -n $n $fastq $quids > >(bgzip > $out/$(basename $fastq))
