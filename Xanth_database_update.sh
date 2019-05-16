#!/bin/bash

for file in ~/Documents/Xanthomonas_blastdb/genomes_to_add/*; do
f=$(sed 's/\(.*\).fasta/\1/' <<< "$(basename $file)")
sed "s/\(>[0-9]*\).*/>\1_$f/" $file > ~/Documents/Xanthomonas_blastdb/genomes_in_db/${f}.dbready.fasta
done

cat ~/Documents/Xanthomonas_blastdb/genomes_in_db/*fasta > ~/Documents/Xanthomonas_blastdb/genomes_in_db/genomes.cat.fasta
mv ~/Documents/Xanthomonas_blastdb/XanthDB* database_backup
makeblastdb -in ~/Documents/Xanthomonas_blastdb/genomes.cat.fasta -dbtype nucl -out ~/Documents/Xanthomonas_blastdb/XanthDB
rm ~/Documents/Xanthomonas_blastdb/genomes_in_db/genomes.cat.fasta
