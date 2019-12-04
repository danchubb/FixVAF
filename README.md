# FixVAF
Code to remove bias from Isaac aligned data by clipping all reads for variant positions by 5 bases and producing a modified vcf file. 

Requires python 3 with psam installed. 

To run:



vcf_file=sys.argv[1]
bam=sys.argv[2]
myfasta=pysam.FastaFile(sys.argv[3])
