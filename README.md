# FixVAF
Code to remove bias from Isaac aligned data by clipping all reads for variant positions by 5 bases and producing a modified vcf file. 

Requires python 3 with psam installed. 

To run:

python FixVaf.py [vcf file] [bam file] [fasta file]

bam file is the file used to call the vcf. 
