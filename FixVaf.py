import pysam
import sys
import numpy
import os
#import matplotlib.pyplot as plt

myfasta=pysam.FastaFile('/scratch/DGE/MOPOPGEN/studies/references/b38/ilumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa')


chromlist=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

chrom_lengths={}

for i in chromlist:
    chrom_lengths[i]=myfasta.get_reference_length(i)

def bamquery(myfasta,bamfile,chr,pos,ref,alt):
    bam=pysam.AlignmentFile(bamfile,"rb")
    b0pos=pos-1
    output={}
    tiers={}
    #var=['ref','alt','other']
    var=['A','T','C','G']
    #fields=['pos','cig','mq','bq']
    fields=['pos','cig','mq','bq','paired']
    tier=['t1','t2','t3']
     
    for v in var:
        tiers[v]={}
        for t in tier:
            tiers[v][t]=0
        for f in fields:
            output[v+'_'+f]=[]
    
    #read position has a huge buffer incase it is clipped at the begining or end. In theory it could be clipped up to 149 bases and then additional indels hsift it further. In reality It shouldn't be that bad but by how much should we shift it? By putting a large window it will capture everything and then ignore everything that does not overlap. I guess it takes longer though so there's probably a better way.
    
    s=b0pos-2
    e=pos+2
    
    if b0pos-2 <0:
        s=0
    if pos+2 > chrom_lengths[chr]:
        e=chrom_lengths[chr]
    readnames={}
    readnames['t1']={}
    readnames['t2']={}
    readnames['t3']={}
    for read in bam.fetch(contig=chr,start=s,end=e):
        cigar=read.cigartuples
        if read.is_duplicate or read.is_unmapped or read.is_secondary or read.is_qcfail or not cigar or read.mapping_quality==0:
            continue
        skip=0
        #skip reads with indels
        #for c in cigar:
        #    if c[0]==2 or c[0]==1:
        #        skip=1
        #if skip:
        #    pass
            #continue
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
                        #placeholder if to keep indentation for when it's ready to be fixed
                        base=read.seq[read_pos]
                        #cigar_val=c[0]
                        bq=read.query_qualities[read_pos]
                        #placeholder for max of trimmed, only of value where there is trimming
                        #max_bq_trimmed=60
                        if(bq>0):
                            pair=0
                            if read.is_proper_pair:
                                pair=1
                            position_in_read=read_pos
                            mq=read.mapping_quality
                            tiers.setdefault(base,{})
                            if cigar_val==0:
                                if mq >= 40 and pair:
                                    if read.query_name in readnames['t1'] and readnames['t1'][read.query_name]==base:
                                        pass
                                    else:
                                        readnames['t1'][read.query_name]=base
                                        tiers[base]['t1']+=1
                                        print(read.query_name,base)
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
                            output[base+'_pos'].append(position_in_read)
                            output[base+'_cig'].append(cigar_val)
                            output[base+'_mq'].append(mq)
                            output[base+'_bq'].append(bq)
                            output[base+'_paired'].append(pair)                                
    return tiers

#sample=sys.argv[1]
vcf_file=sys.argv[1]
bam=sys.argv[2]
#strelka_bam=sys.argv[3]
#bwa_bam=sys.argv[4]


def read_somatic_vcf(vcf_file,bam):
    vcf_out=vcf_file[:-7]+".VAF.vcf"
    OUT=open(vcf_out,'w')
    vcf=pysam.VariantFile(vcf_file,'r')
    sample_list=list(vcf.header.samples)
    #vcf=pysam.VariantFile('/data/rds/shared/MYELOMA/Ilumina_WGS/38397484/Variations/LP2100049-DNA_A01.vcf.gz','r')
    
    
    counts=0
    
    vcf.header.add_line('##INFO=<ID=t1BamAC,Number=1,Type=String,Description="Allele counts of variant from Bamfile after all reads softclipped by at least 5BP from ends. MQ>40, Reads must be in proper pair, BQ>0. Order is A,C,G,T">')
    vcf.header.add_line('##INFO=<ID=t2BamAC,Number=1,Type=String,Description="Allele counts of variant from Bamfile after all reads softclipped by at least 5BP from ends. MQ>5, BQ>0. Order is A,C,G,T">')
    vcf.header.add_line('##INFO=<ID=t3BamAC,Number=1,Type=String,Description="Allele counts of variant from Bamfile after all reads softclipped by at least 5BP from ends. MQ>0, BQ>0. Order is A,C,G,T">')
    #print(str(vcf.header).rstrip())
    OUT.write(str(vcf.header))
    #182791
    for v in vcf:
        filters=v.filter.keys()
        sample=v.samples[1]
        #GT=v.samples[sample_list[0]]['GT']
        #print(v)
        #if filters[0]=="PASS" and v.chrom in chromlist and len(v.alts)==1 and len(v.ref)==1 and len(v.alts[0])==1 and len(GT)==2 and sum(GT)==1:
        if filters[0]=='PASS' and v.chrom in chromlist and v.alts and len(v.alts)==1 and 'TIR' not in v.format :
            #(len(v.ref)>1 or len(v.alts[0])>1):
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


read_somatic_vcf(vcf_file,bam)
