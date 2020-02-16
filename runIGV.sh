# files needed in my computer
# /mnt/d/GitHub/Course7Tutor/SD.fa
# /mnt/d/GitHub/Course7Tutor/output.sorted.bam
# Indices for previous
# /mnt/d/GitHub/Course7Tutor/output.sorted.bam.bai

printf "moving to IGV directory\n" &&
cd igv/IGV_2.8.0/ &&
printf "running IGV\n" &&
java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main
