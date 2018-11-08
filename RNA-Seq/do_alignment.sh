#!/bin/bash -e

# do alignment for pair end fastq.gz file,
# and sort, index the result bam file
# 
# aligner:
#   bwa
#       aln
#       mem
#   bowtie2
#   hisat2
#
# fastq should named as _R1.fastq.gz _R2.fastq.gz


## CONFIG

idx_prefix="/media/nanguage/DATA/Genomes/MTB/bwa_index/H37Ra"

reads=$(ls ./*.fastq.gz)

aligner=bwa-aln # bowtie2 / bwa-aln / bwa-mem / hisat2

threads=4

pbs=0 # `1` for use pbs system

## END CONFIG


function align{
    id=$1

    echo $id
    echo do alignment with $aligner

    R1="$id"_R1.fastq
    R2="$id"_R2.fastq

    if [ "$aligner" == "bwa-aln" ]; then
        echo "decompressing ..."
        zcat $R1.gz > $R1 &
        zcat $R2.gz > $R2 &
        wait

        echo "alignment ..."
        bwa aln $idx_prefix -t $threads $R1 > "$id"_R1.sai
        bwa aln $idx_prefix -t $threads $R2 > "$id"_R2.sai

        echo "sai to sam ..."
        bwa sampe $idx_prefix "$id"_R1.sai "$id"_R1.sai $R1 $R2 > $id.sam
        rm "$id"_R1.sai
        rm "$id"_R2.sai
    elif [ "$aligner" == "bwa-mem" ]; then
        echo "decompressing ..."
        zcat $R1.gz > $R1 &
        zcat $R2.gz > $R2 &
        wait

        echo "alignment ..."
        bwa mem $idx_prefix $R1 $R1 -t $threads > $id.sam
    elif [ "$aligner" == "bowtie2" ]; then
        bowtie2 -x $idx_prefix -1 $R1 -2 $R2 -S $id.sam -p $threads
    elif [ "$aligner" == "hisat2" ]; then
        hisat2 -x $idx_prefix -1 $R1 -2 $R2 -S $id.sam -p $threads
    fi

    echo "samtools sort"
    samtools view -bh $id.sam > $id.bam
    samtools sort -@ $threads $id.bam -o $id.sorted.bam
    rm $id.sam $id.bam
    samtools index -@ $threads $id.sorted.bam

    echo "done"
    rm $R1 $R2
}


export -f align


# do all align task
for id in `echo "$reads" | sed 's/.fastq.gz//g' | sed 's/_R1//' | sed 's/_R2//' | sort -u`;
do
    if [ $pbs -eq 1 ]; then
        echo "align $id" | qsub -V -l nodes=1:ppn=$threads -N align_"$id" -d $PWD
    else
        align $id
    fi
done

