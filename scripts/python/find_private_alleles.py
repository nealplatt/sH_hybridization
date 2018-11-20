import sys

#open necessary files (read in from the cmd line)
#vcf_in=open(sys.argv[1], 'r')
#vcf_out=open(sys.argv[2], 'w')
vcf_in.close()
snps_out.close()

vcf_in=open("auto_maf.vcf", 'r')
snps_out=open("snps.out", 'w')

informative_snps={}

string='\t'.join(["chrom", "pos", "id", "bovis_0", "tz_0", "out_0", "bovis_1", "tz_1", "out_1"])
snps_out.write(string + '\n')

for vcf_entry in vcf_in:
    
    if not vcf_entry.startswith("#"):
        vcf_entry=vcf_entry.rstrip()
        #its an entry and the alt/ref needs to be swapped and 0/1s tr'd
        (chrom, pos, id, ref, alt, qual, filter, info, format)=vcf_entry.split("\t")[0:9]
        gt=vcf_entry.split("\t")[9:]
        
        #get the alleles for each "pop"
        b=[i.split(':')[0] for i in gt][1]
        o=[i.split(':')[0] for i in gt][2:8]
        n=[i.split(':')[0] for i in gt][9:54]
        t=[i.split(':')[0] for i in gt][55:101]
        
        oa=list(set([i.split(':')[0] for i in gt][2:8]))
        na=list(set([i.split(':')[0] for i in gt][9:54]))
        ta=list(set([i.split(':')[0] for i in gt][55:101]))
        
        if b is not './.':
           
            tas={'1':0, '0':0}
            for i in t:
                tas[i.split('/')[0]]=1
                tas[i.split('/')[1]]=1
            
            bas={'1':0, '0':0}
            for i in b.split('/'):
                bas[i]=1
            
            oas={'1':0, '0':0}
            for i in o:
                oas[i.split('/')[0]]=1
                oas[i.split('/')[1]]=1
            
            nas={'1':0, '0':0}
            for i in n:
                nas[i.split('/')[0]]=1
                nas[i.split('/')[1]]=1    
            
            
            #identify ancestry informative snps            
            if bas['0'] == 1 and tas['0'] == 0 and oas['0'] == 0:
                #anc_inf="bovis-0"
                #snps_out.write(string + '\t' + anc_inf + '\n')
                informative_snps[id]="b;0"
            elif bas['1'] == 1 and tas['1'] == 0 and oas['1'] == 0:
                #anc_inf="bovis-1"
                #snps_out.write(string + '\t' + anc_inf + '\n')
                informative_snps[id]="b;1" 
            elif tas['0'] == 1 and bas['0'] == 0 and oas['0'] == 0:
                #anc_inf="tz-0"
                #snps_out.write(string + '\t' + anc_inf + '\n')
                informative_snps[id]="t;0"    
            elif tas['1'] == 1 and bas['1'] == 0 and oas['1'] == 0:
                #anc_inf="tz-1"
                #snps_out.write(string + '\t' + anc_inf + '\n')
                informative_snps[id]="t;1"  
            else:
                #informative_snps[id]="a"
                pass
            
            string='\t'.join([chrom, pos, id, str(bas['0']), str(tas['0']), str(oas['0']), str(bas['1']), str(tas['1']), str(oas['1'])])
            snps_out.write(string + '\n')
            
            b=[]
            o=[]
            n=[]
            t=[]
            anc_inf="null"
            
            #now use this info to modify the vcf file to contain snp ancestry

vcf_in.close()
snps_out.close()


#now that we have a list of informative snps. read the phased file and
#H1 TTTTBBTBTBBTBBTBBBBBBTTT
#h2 TTTTTTTTTTTTTTTTTTTBBBTT

phased_vcf_in=open("auto_beagle_maf05.vcf", 'r')

positions={}
chromosomes={}
haplotypes={}        
snp_order=[]

for phased_vcf_entry in phased_vcf_in:
    if phased_vcf_entry.startswith("#") and not phased_vcf_entry.startswith("##"):
        header_line=phased_vcf_entry.rstrip()
        sample_ids=header_line.split("\t")[9:]
    
    
    if not phased_vcf_entry.startswith("#"):
        phased_vcf_entry=phased_vcf_entry.rstrip()
        #its an entry and the alt/ref needs to be swapped and 0/1s tr'd
        (chrom, pos, id, ref, alt, qual, filter, info, format)=phased_vcf_entry.split("\t")[0:9]
        
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


