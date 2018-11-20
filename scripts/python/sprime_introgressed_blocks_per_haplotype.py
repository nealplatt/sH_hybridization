import sys

#open necessary files (read in from the cmd line)
#vcf_in=open(sys.argv[1], 'r')
#vcf_out=open(sys.argv[2], 'w')
vcf_in.close()
sprime_scores_in.close()
initial_table_out.close()

vcf_in=open("sprime_snps_maf05.vcf", 'r')
sprime_scores_in=open("niger_tz-out_sprime_2018-06-12.score", 'r')
initial_table_out=open("initial_table.out", 'w')
#snps_out=open("snps.out", 'w')

first_line = True
introgressed_alleles ={}
introgressed_loci = []

for entry in sprime_scores_in:
    if first_line:
        first_line = False
        continue    
    entry=entry.rstrip()
    (chrom, pos, id, ref, alt, segment, allele, score)=entry.split("\t")
        
    introgressed_loci.append(id)
    introgressed_alleles[id]=allele
    

modified_genotypes=[] 
matrix=[]   
num_intro = 0
for vcf_entry in vcf_in:
    
    if not vcf_entry.startswith("#"):
        vcf_entry=vcf_entry.rstrip()
        
        #its an entry and the alt/ref needs to be swapped and 0/1s tr'd
        (chrom, pos, id, ref, alt, qual, filter, info, format)=vcf_entry.split("\t")[0:9]
        
        genotypes=[i.split(':')[0] for i in vcf_entry.split("\t")[9:]]
        
        if id in introgressed_loci:
            introgressed_allele=introgressed_alleles[id]        
            
            if introgressed_allele == "0":
                genotypes = [gt.replace('0', 'I') for gt in genotypes]
                genotypes = [gt.replace('1', 'N') for gt in genotypes]
            
            elif introgressed_allele == "1":
                genotypes = [gt.replace('1', 'I') for gt in genotypes]
                genotypes = [gt.replace('0', 'N') for gt in genotypes]
             
            h_a=([i.split('|', 1)[0] for i in genotypes])
            h_b=([i.split('|', 1)[1] for i in genotypes]) 
            
            entry='\t'.join([chrom, pos, id, ref, alt, introgressed_allele]+genotypes)
            matrix.append([chrom, pos, id, ref, alt, introgressed_allele] + h_a + h_b)
        else:
            
            genotypes = [gt.replace('0', 'A') for gt in genotypes]
            genotypes = [gt.replace('1', 'A') for gt in genotypes]
            
            h_a=([i.split('|', 1)[0] for i in genotypes])
            h_b=([i.split('|', 1)[1] for i in genotypes])
            
            entry='\t'.join([chrom, pos, id, ref, alt, "A"]+genotypes)
            matrix.append([chrom, pos, id, ref, alt, "A"] + h_a + h_b)
        
        initial_table_out.write(entry + '\n')
        modified_genotypes.append(entry)



vcf_in.close()
initial_table_out.close()

#get sample ids
vcf_in=open("sprime_snps_maf05.vcf", 'r')
for vcf_entry in vcf_in:
    
    if vcf_entry.startswith("#CHROM"):
        vcf_entry=vcf_entry.rstrip()
        sample_ids=[i.split(':')[0] for i in vcf_entry.split("\t")[9:]]

vcf_in.close()

haplotype_ids=[]
for sample in sample_ids:
    haplotype_ids.append(sample + "_A")
    haplotype_ids.append(sample + "_B")


################################################################################

#create strings for each sample/haplotype
for locus in matrix
    haplotypes = matrix[1][6:] 
    for haplotype in list(range(0:len(haplotypes))):
        sequence[haplotype_id]=haplotypes[haplotype] 
           


#convert matrix to numpy array
#iterate through each of the "haplotype columns"
#convert to string

import numpy as np
np_matrix=np.array(matrix)

np_matrix

for i in list(range(6,197)):
    x=''.join(np_matrix[:,i])
    print(x)










string="NNNNNNNIIIIINNNNINNNNNNIIIIIIIIIIIIINNINIINI"
[m.span() for m in re.finditer("II+", string)]     


positions={}
chromosomes={}
haplotypes={}        
snp_order=[]

matrix= []
h_a=[]
h_b=[]

for entry in modified_genotypes:
    (chrom, pos, id, ref, alt, introgressed_allele)=entry.split("\t")[0:6]
    genotypes=[i.split(':')[0] for i in entry.split("\t")[6:]]
    h_a=([i.split('|', 1)[0] for i in genotypes])
    h_b=([i.split('|', 1)[1] for i in genotypes])
    
    matrix.append([chrom, pos, id, ref, alt, introgressed_allele] + h_a + h_b)
    
first_line = True
for sample in list(range(0, len(sample_ids))):
    if first_line == True:
        for haplotype in h_a:
            a_haplotypes[sample_ids[sample]]=[h_a[sample]]
            b_haplotypes[sample_ids[sample]]=[h_b[sample]]
            print("here")
        first_line = False


    else:
        for haplotype in a_haplotypes:
            a_haplotypes{sample_id[sample]}.append[h_a[sample]]
            b_haplotypes{sample_id[sample]}.append[h_a[sample]]
      



        if id in informative_snps:
            gt=phased_vcf_entry.split("\t")[9:]
            
            gt_a=([i.split('|', 1)[0] for i in gt])
            gt_b=([i.split('|', 1)[1] for i in gt])
            
            (anc, allele)=informative_snps[id].split(';')
            if anc == 'b' and allele == '0':
                gt_a=['b' if i == '0' else "a" for i in gt_a] 
                gt_b=['b' if i == '0' else "a" for i in gt_b]
            elif anc == 'b' and allele == '1':
                gt_a=['b' if i == '1' else "a" for i in gt_a]
                gt_b=['b' if i == '1' else "a" for i in gt_b] 
            elif anc == 't' and allele == '0':
                gt_a=['t' if i == '0' else "a" for i in gt_a] 
                gt_b=['t' if i == '0' else "a" for i in gt_b]
            elif anc == 't' and allele == '1':
                gt_a=['t' if i == '1' else "a" for i in gt_a] 
                gt_b=['t' if i == '1' else "a" for i in gt_b]
            else:
                pass
            
            positions[id]=pos
            chromosomes[id]=chrom
            haplotypes[id]=gt_a + gt_b
            
            snp_order.append(id)

phased_vcf_in.close()

anc_out=open("anc.tsv", 'w')
for i in range(len(snp_order)):
    snp=snp_order[i]
    
    t_haplotypes='\t'.join(haplotypes[snp])
    
    anc_out.write(chromosomes[snp] + '\t' + positions[snp] +  '\t' + snp + '\t' + t_haplotypes + '\n')

anc_out.close()

#THERE ARE MANY FEWERE SNPS IN OUTPUT FILE THAN EXPECTED FROM CHECK INTO THIS
# EX
#len(informative_snps)
#31166
#len(snp_order)
#17457

#why are there so many fewer....part of it has to do with fewer snps in the beagle file

#for each individual find consequtive B alleles then print out start and stop


#get names for blocks file
haplotype_ids=[]
for a in sample_ids:
    haplotype_ids.append(a + '_A')

for b in sample_ids:
    haplotype_ids.append(b + '_B')


blocks_out=open("blocks.tsv", 'w')
#write header to blocks file
blocks_out.write("indiv_hap" + '\t' "chrom" + '\t' + "Start_bp" + '\t' + "end_bp" + '\t' + "length_bp" + '\t' + "start_snp" + '\t' + "end_snp" + '\t' + "snps_in_block" + '\n')

#cycle through each individual
for i in range(len(matrix[1])):
    h=[]
    locus=0
    chain_start=0
    chain_end=0
    chain_length=0   
    
    hap_id=haplotype_ids[i]
    #build an array of snps for that individual    
    for j in range(len(matrix)):                                                                                                                
        h.append(matrix[j][i])
    
    print(str(len(h)))
    #iterate over the array looking for blocks of "b[b|a]*b"
    while locus < len(h):
        next=1
        if h[locus] is 'b':
        #start the chain
            #print("Starting new chain at: " + str(locus) + ": " + h[locus])
            chain_start=locus
            chain_end=locus
            while locus+next < len(h) and h[locus+next] is not 't' and  chromosomes[snp_order[locus]] == chromosomes[snp_order[locus+next]]:
                if h[locus+next] is 'b':
                    chain_end=locus+next
                #extend the chain one more
                #print("...ext: " + str(locus+next) + ": " + h[locus+next])
                next=next+1
            #print("...term: " + str(locus+next) + ": " + h[locus+next])
            if chain_start is not chain_end:
                chromosome = chromosomes[snp_order[chain_start]]
                chain_start_snp = snp_order[chain_start]
                chain_end_snp = snp_order[chain_end]
                chain_start_bp = positions[snp_order[chain_start]]
                chain_end_bp = positions[snp_order[chain_end]]
                chain_length = str(int(chain_end)-int(chain_start)+1)
                chain_length_bp = str(int(chain_end_bp) - int(chain_start_bp) + 1) 
                
                out='\t'.join([hap_id, chromosome, chain_start_bp, chain_end_bp, chain_length_bp, chain_start_snp, chain_end_snp, chain_length])
                #print(out)
                blocks_out.write(out + '\n')
                            
                chain_start=0
                chain_end=0
                chain_length=0
                
        locus=locus+next

blocks_out.close()

bovis_snps_out=open("bovis_snps.txt", 'w')
#get bovis autapomorphic loci
for snp in informative_snps:
    if informative_snps[snp].split(';')[0] is 'b':
         bovis_snps_out.write(snp +'\t' + chromosomes[snp] +'\t' + positions[snp] +'\n')

bovis_snps_out.close()

#
vcftools \
    --vcf auto_beagle_maf05.vcf \
    --snps bovis_snps.list \
    --keep niger.txt \
    --recode \
    --stdout \
    >ne_autapomoprhic.vcf

plink \
    --vcf ne_autapomoprhic.vcf \
    --out ne_autapomoprhic_ld \
    --r dprime bin\
    --allow-extra-chr

vcftools \
    --vcf ne_autapomoprhic.vcf \
    --hap-r2 \
    --stdout \
    >ne_autapomoprhic.ld


