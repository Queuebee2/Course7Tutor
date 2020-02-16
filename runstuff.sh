#!/usr/bin/env bash

printf "be aware that to use this script, Xming or x11 or something must be running for IGV to start at the end\n"
sleep 3
printf "Creating indexes from reference genome ("SD.fa") with bowtie2-build\n" &&
bowtie2-build "SD.fa" SD &&
printf "Indexes successfully made\n" &&
printf "Creating a sequence alignment map with bowtie2  \n" &&
bowtie2 --phred64 --very-fast --no-unal -x SD -1 SD_3_1_400bp.fq -2 SD_3_2_400bp.fq -S Docc.sam &&
printf "successfully ran bowti2\n" &&
printf "creating view with samtools\n" &&
samtools view -b -o Docc.sort.bam -S Docc.sam &&
printf "successfully created view with samtools\n" &&
samtools sort Docc.sort.bam > Docc.sorted.bam &&
printf "sorted with samtools\n"  &&
samtools index Docc.sorted.bam &&
printf "successfully created BAI (indexing) with samtools\n" &&
samtools depth Docc.sorted.bam > Docc.sorted.depth &&
printf "successfully did depth thingy with samtools\n" &&
printf "moving to IGV directory\n" &&
cd igv/IGV_2.8.0/ &&
printf "running IGV\n" &&
java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main
