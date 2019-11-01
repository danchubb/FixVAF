import pysam
import sys
import numpy
import os
#import matplotlib.pyplot as plt
#Also requires tabix/bgzip to be installed. 

    
def input_output_same(input,output):
    intype='cat'
    outtype='cat'
    
    if input[-2:]=='gz':
        intype='z'+intype
    if output[-2:]=='gz':
        outtype='z'+outtype
    if os.path.exists(input) and os.path.exists(output):
        invars=int(os.open(intype+' '+input+' |grep -v "#" |wc').read().split()[0])
        outvars=int(os.open(outtype+' '+output+' |grep -v "#" |wc').read().split()[0])
        if invars == outvars:
            return 0
        else:
            print("variants in "+input+" ("+str(invars)+") not equal to "+output+"("+str(outvars)+")",file=sys.stderr)
            return 1
    else:
        print("Check input and output files",input,output,file=sys.stderr)
        return 1

#'/scratch/DGE/MOPOPGEN/studies/references/b38/ilumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa')
vcf_file=sys.argv[1]
bam=sys.argv[2]
myfasta=pysam.FastaFile(sys.argv[3])

chromlist=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
chrom_lengths={}
for i in chromlist:
    chrom_lengths[i]=myfasta.get_reference_length(i)

def bamquery(myfasta,bamfile,chr,pos,ref,alt):
    bam=pysam.AlignmentFile(bamfile,"rb")
    b0pos=pos-1
    tiers={}
    var=['A','T','C','G']
    fields=['pos','cig','mq','bq','paired']
    tier=['t1','t2','t3']
    readnames={} 
    for v in var:
        tiers[v]={}
        for t in tier:
            tiers[v][t]=0
            readnames[t]={}
    
    s=b0pos-2
    e=pos+2
    
    if b0pos-2 <0:
        s=0
    if pos+2 > chrom_lengths[chr]:
        e=chrom_lengths[chr]
    
    for read in bam.fetch(contig=chr,start=s,end=e):
        cigar=read.cigartuples
        #ignore unmapped, duplicate, non primary badly mapped or cigar is malformed 
        if read.is_duplicate or read.is_unmapped or read.is_secondary or read.is_qcfail or not cigar or read.mapping_quality==0:
            continue
        #where there is a hardclip, begining, end or both? Need to know in the are cases where there is both soft and hard clipping happening
        H_begin=0
        H_end=0
        hard_count=0
        for c in cigar:
            hard_count+=1
            if c[0]==5:
                if hard_count==1:
                    H_begin=1
                elif hard_count==len(cigar):
                    H_end=1
        cig_hash={}
        c_pos=read.pos
        
        read_pos=0
            
        seq_map={}
        
        cig_place=0
        base=''
        for c in cigar:
            #if we haven't found the base yet
            if not base:
                cig_place+=1
                if c[0] == 5:
                    continue
                elif(c[0] == 1):
                    read_pos+=c[1]
                elif c[0] == 2:
                    c_pos+=c[1]
                elif c[0]== 4 or c[0]==0:
                    #if soft clip is of 5' then adjust alignment otherwise don't adjust because it's at the end and is correctly placed. The position of the end will be size of cigar - potential following hardclip
                    if c[0] == 4 and cig_place!=len(cigar)-H_end:
                        c_pos+=-c[1]
                        read.pos+=-c[1]
                    if b0pos>c_pos+c[1]-1:
                        read_pos+=c[1]
                        c_pos+=c[1]
                    elif b0pos>=c_pos and b0pos<=c_pos+c[1]-1:
                        diff=b0pos-c_pos
                        c_pos+=diff
                        read_pos+=diff
                        read_size=len(read.seq)
                        cigar_val=c[0]
                        if read_pos<=4 or read_pos>=(read_size-5):
                            cigar_val=4
                        base=read.seq[read_pos]
                        bq=read.query_qualities[read_pos]
                        if(bq>0):
                            pair=0
                            if read.is_proper_pair:
                                pair=1
                            position_in_read=read_pos
                            mq=read.mapping_quality
                            tiers.setdefault(base,{})
                            if cigar_val==0:
                                
                                #define the three tiers. As with allelecounter, if the same position is overlapped by both ends of a read then count it once if it is the same base and twice if it is a different one. Alternative would be to ignore entirely but wither way there are edge cases
                                
                                if mq >= 40 and pair:
                                    if read.query_name in readnames['t1'] and readnames['t1'][read.query_name]==base:
                                        pass
                                    else:
                                        readnames['t1'][read.query_name]=base
                                        tiers[base]['t1']+=1
                                if mq >= 5:
                                    if read.query_name in readnames['t2'] and readnames['t2'][read.query_name]==base:
                                        pass
                                    else:
                                        readnames['t2'][read.query_name]=base
                                        tiers[base]['t2']+=1
                                if mq > 0:
                                    if read.query_name in readnames['t3'] and readnames['t3'][read.query_name]==base:
                                        pass
                                    else:                                
                                        tiers[base]['t3']+=1
                                        readnames['t3'][read.query_name]=base
                        
    return tiers

def read_somatic_vcf(vcf_file,bam):
    vcf_out=vcf_file[:-7]+".VAF.vcf"
    OUT=open(vcf_out,'w')
    vcf=pysam.VariantFile(vcf_file,'r')
    sample_list=list(vcf.header.samples)
   
    counts=0
    
    vcf.header.add_line('##INFO=<ID=t1BamAC,Number=4,Type=Integer,Description="Allele counts of variant from Bamfile after all reads softclipped by at least 5BP from ends. MQ>40, Reads must be in proper pair, BQ>0. Order is A,C,G,T">')
    vcf.header.add_line('##INFO=<ID=t2BamAC,Number=4,Type=Integer,Description="Allele counts of variant from Bamfile after all reads softclipped by at least 5BP from ends. MQ>5, BQ>0. Order is A,C,G,T">')
    vcf.header.add_line('##INFO=<ID=t3BamAC,Number=4,Type=Integer,Description="Allele counts of variant from Bamfile after all reads softclipped by at least 5BP from ends. MQ>0, BQ>0. Order is A,C,G,T">')
    OUT.write(str(vcf.header))
    for v in vcf:
        filters=v.filter.keys()
        sample=v.samples[1]
        #if filters[0]=="PASS" and v.chrom in chromlist and len(v.alts)==1 and len(v.ref)==1 and len(v.alts[0])==1 and len(GT)==2 and sum(GT)==1:
        if v.chrom in chromlist and v.alts and len(v.alts)==1 and len(v.ref)==1 and len(v.alts[0])==1 and 'TIR' not in v.format :
            tiers=bamquery(myfasta,bam,v.chrom,v.pos,v.ref,v.alts[0])
            v.info['t1BamAC']=str(tiers['A']['t1'])+","+str(tiers['C']['t1'])+","+str(tiers['G']['t1'])+","+str(tiers['T']['t1'])
            v.info['t2BamAC']=str(tiers['A']['t2'])+","+str(tiers['C']['t2'])+","+str(tiers['G']['t2'])+","+str(tiers['T']['t2'])
            v.info['t3BamAC']=str(tiers['A']['t3'])+","+str(tiers['C']['t3'])+","+str(tiers['G']['t3'])+","+str(tiers['T']['t3'])
            
        vline=str(v)
        vline=vline.rstrip()
        OUT.write(vline+"\n")
    OUT.close()
    os.system("bgzip -f "+vcf_out)
    os.system("tabix -f -p vcf "+vcf_out+".gz")
    same=input_output_same(vcf_file,vcf_out+".gz")
    return same


error=read_somatic_vcf(vcf_file,bam)

if error:
    print("ERROR, number of variants is different or files don't exist", file=sys.stderr)
    sys.exit(1)
else:
    print("Input and output exists and number of vars is the same")
