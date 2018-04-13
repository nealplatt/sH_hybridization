#VCF TO DIPLOID FASTA
#python vcf_to_diploid_fasta.py <input vcf>
import sys

#to do
# check on the sequences.fromsequence keys code
# exclude things that aren't biallelic?
# print to fasta

#have list of ambig codes

homozygous_reference_genotype="0/0"
homozygous_alternate_genotype="1/1"
heterozygous_genotype="0/1"
missing_genotype="./."


iupac_code ={"AC":"M", "AG":"R", "AT":"W", "CG":"S", "CT":"Y", "GT":"K", "CA":"M", "GA":"R", "TA":"W", "GC":"S", "TC":"Y", "TG":"K"}

sequences={}
#open necessary files (read in from the cmd line)
#vcf_file_in=open(sys.argv[1], 'r')
#fasta_file_out=open(sys.argv[2], 'w')

vcf_file_in=open("test.vcf", 'r')

for vcf_entry in vcf_file_in:
    
    vcf_entry=vcf_entry.rstrip()
    
    if vcf_entry.startswith("#CHROM"):
        #its the header line
     
        #create list of samples
        samples = vcf_entry.split("\t")[9:]
        
        #CHECK ON THIS I DON"T UNDERSTAND IT       
        sequences=sequences.fromkeys(samples, '')
        
        #print("if")
    
    elif vcf_entry.startswith("##"):
        a=1
        #print("elif")
    else:
        ref_allele = vcf_entry.split("\t")[3]
        alt_allele = vcf_entry.split("\t")[4]
        genotypes = vcf_entry.split("\t")[9:]
              
        print("else")
        
        for i in range(1, len(genotypes)+1):
            print("for loop")
            sample = samples[i]
            genotype = genotypes[i].split(":")[0]
            
            if genotype == homozygous_reference_genotype:
                nuc = ref_allele
            elif genotype == homozygous_alternate_genotype:
                nuc = alt_allele
            elif genotype == heterozygous_genotype:
                nuc = iupac_code[ref_allele + alt_allele]
            elif genotype == missing_genotype:
                nuc = "?"
            else:
                #print("Unexpected genotype: ", genotype, "...exiting")
                sys.exit
                        
            sequences[sample] = sequences[sample]+nuc


#now pring everything in the sequences dictionary
                



print as fasta
