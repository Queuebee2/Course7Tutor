#!/usr/bin/env bash

#assuming ORFfinder, Repeatscout and augustus are installed

# install ORFfinder
####################
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.g
# gunzip ORFfinder.gz
####################
# thats it

# todo: install repeatscout, augustus


GENOME_FILE="SD.fa"
ORF_FINDER_OUT="found_orfs"
MIN_LENGTH=150
LMER_OUTPUT="output_lmer.freq"
RS_OUTPUT="output_repeatscout.fas"
RS_PARAMETERS=(-minthresh 1 -minimprovement 1 -stopafter 1000)
AUG_SPECIES="debaryomyces_hansenii"
AUG_OUTPUT="augustus_output.gff"

printf "running ORFfinder on $GENOME_FILE with minimum length $MIN_LENGTH\n"
./ORFfinder -in "$GENOME_FILE" -ml "$MIN_LENGTH" -out "$ORF_FINDER_OUT" -outfmt 1

printf "Creating l-mer table with repeatscout build_lmer \ninput sequence=$GENOME_FILE \noutput_file:=$LMER_OUTPUT\n"
./RepeatScout-1/build_lmer_table -sequence "$GENOME_FILE" -freq "$LMER_OUTPUT"

printf "Finding repeats with repeatscout on $GENOME_FILE \nl-mer_frequencies=$LMER_OUTPUT \noutput_file:$RS_OUTPUT\n"
printf "extra parameters:\n"
printf '%s %s \n' "${RS_PARAMETERS[@]}"
./RepeatScout-1/RepeatScout -sequence SD.fa -output output_repeatscout -freq output_lmer.freq "${RS_PARAMETERS[@]}"

printf "Predicting genes with augustus on $GENOME_FILE \nspecies=$AUG_SPECIES \noutput_file:$AUG_OUTPUT\n"
augustus "--species=$AUG_SPECIES" "$GENOME_FILE" > "$AUG_OUTPUT" && echo "done"
