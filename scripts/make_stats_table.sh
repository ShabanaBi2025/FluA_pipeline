#!/bin/bash

# Header
echo -e "Sample\tRecords\tSNPs\tIndels\tTs\tTv\tTs/Tv\tSingletons\tMost Frequent QUAL Bin"

# Loop through all .vcf.stats files
for stats in *.vcf.stats; do
    sample=$(basename "$stats" .vcf.stats)

    records=$(awk '/number of records:/ {print $NF}' "$stats")
    snps=$(awk '/number of SNPs:/ {print $NF}' "$stats")
    indels=$(awk '/number of indels:/ {print $NF}' "$stats")
    ts=$(awk '/^TSTV/ {print $3}' "$stats")
    tv=$(awk '/^TSTV/ {print $4}' "$stats")
    ratio=$(awk '/^TSTV/ {print $5}' "$stats")
    singletons=$(awk '/^SiS/ {print $4}' "$stats")
    qual_bin=$(awk '/^QUAL/ {print $4 "\t" $5}' "$stats" | sort -k2 -nr | head -1 | awk '{print $1 " (" $2 " SNPs)"}')

    echo -e "$sample\t$records\t$snps\t$indels\t$ts\t$tv\t$ratio\t$singletons\t$qual_bin"
done
