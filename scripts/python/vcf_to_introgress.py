
!#/master/nplatt/miniconda3/envs/snp_calling/bin/python
# vcf_to_introgress.py <in_vcf> <introgress_geno_matrix> <introgress_loci>

#this program takes in a vcf file and converts it into X csv files for
# import into the R introgress package.  

#IT IS VERY IMPORTANT TO NOTE THAT BOVIS IN THIS CASE IS THE SECOND SAMPLE
# IN THE VCF

# also the matrix will have a B for the bovis allele and an H for the haem
#  allele
import sys

#vcf_file_in=open("sys.argv[1]", 'r')
#matrix_file_out=open("sys.argv[2]", 'w')
#loci_file_out=open("sys.argv[3]", 'w')

vcf_file_in=open("bovis_zanzibar_fixed_snps_autosomal_LD.vcf", 'r')
matrix_file_out=open("matrix", 'w')
loci_file_out=open("loci", 'w')


bov_11_intab="10"
bov_11_outtab="BH"
bov_00_intab="01"
bov_00_outtab="BH"


#alter the chr and start position with the 
#  information from the bed file (stored as a dict)
for vcf_entry in vcf_file_in:
    vcf_entry=vcf_entry.rstrip()
    if not vcf_entry.startswith("#"):
        locus_info = vcf_entry.split("\t")
        gt_info = locus_info[9:]
        locus_name = locus_info[2]
        
        ref_allele=locus_info[3]
        alt_allele=locus_info[4]
        #cycle through for each sample gt and clean it up        
        
        gt=[]
        for indiv_gt in gt_info:
            gt.append(indiv_gt.split(":")[0])
            
        cleaned_gt=",".join(gt)
        
        bovis_gt=cleaned_gt.split(",")[1]
        if bovis_gt == '0/0':
            
            trantab=cleaned_gt.maketrans(bov_00_intab, bov_00_outtab)
            cleaned_gt=cleaned_gt.translate(trantab)
            
            #cleaned_gt=cleaned_gt.replace('B', 'P2')
            #cleaned_gt=cleaned_gt.replace('H', 'P1')
            #change all of the 0s to P2
        elif bovis_gt == '1/1':
            
            trantab=cleaned_gt.maketrans(bov_11_intab, bov_11_outtab)
            cleaned_gt=cleaned_gt.translate(trantab)
            
            #cleaned_gt=cleaned_gt.replace('H', 'P1')
            #cleaned_gt=cleaned_gt.replace('B', 'P2')
        else:
            print("something bad happened:", bovis_gt)
            #continue
                
        cleaned_gt=cleaned_gt.replace('.', 'NA')
        cleaned_gt=cleaned_gt.replace("B/H", "H/B")        
        
        #print gt matrix        
        matrix_file_out.write("%s\n" % cleaned_gt)
        #print locus matrix
        loci_file_out.write("%s,%s\n" % (locus_name, "h"))

vcf_file_in.close()
matrix_file_out.close()
loci_file_out.close()

#added iniv and pop information post script

