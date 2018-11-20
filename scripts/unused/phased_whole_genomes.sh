source activate snp_calling
source /master/nplatt/schisto_hybridization/scripts/set_env.sh


cd $RESULTS_DIR
mkdir invadolysin_gene_family_evol
cd invadolysin_gene_family_evol

#mask genome for unsampled regions

#first convert all bams to bed and merge

SAMPLES=(   "Sh_Dai_044_1"          "Sh_DaiCP_276_1"        "Sh_Kar_001_1"
            "Sh_Kar_37_2"           "Sh_Lata_078_1"         "Sh_PEM_0103_1"
            "Sh_PEM_0104_1"         "Sh_PEM_0130_1"         "Sh_Tiag_272_1"
            "Sh_UNG_0038_1"         "Sh_UNG_0099_1"         "Sh_UNG_0121_1" )


#get sHaem annotation from schisto db

#find invadolysin/leishmanolysin IDs

grep -i leish SchistoDB-40_ShaematobiumEgypt.gff \
    | grep gene \
    | cut -f9 \
    | cut -f1 -d";" \
    | sort -u \
    | sed 's/ID=//' \
    >leishmanolysin_ids.list 




#get invadolysin genes of interest...first need to update old SMPs with new
while read GENE; do
    grep $GENE SchistoDB-40_ShaematobiumEgypt.gff \
        | grep "CDS" >$GENE.gff

    awk '{SUM+=$5-$4+1} END {print SUM/3}' $GENE.gff
done <leishmanolysin_ids.list 

#this gives a list of the cds's for each gene (they all are in frame)
#now check coverage for these regions from the whole genome data.
cat MS*gff | awk '{print $1".1\t"$4"\t"$5}' >schHae_leish_cds.bed

#create masking file for each sample


#for each individual

#calculate coverage
bedtools coverage \
        -a schHae_leish_cds.bed \
        -b ../map_reads/exome/ERR103048_processed.bam \
        -d \
        >ERR103048_leish_cds_cov.bed &

for SAMPLE in "${SAMPLES[@]}"; do
    bedtools coverage \
        -a schHae_leish_cds.bed \
        -b ../map_reads/wgrs/"$SAMPLE"_wgrs_processed.bam \
        -d \
        >"$SAMPLE"_leish_cds_cov.bed
done


INDIVIDUALS=(    "Sh_Dai_044_1"          "Sh_DaiCP_276_1"        "Sh_Kar_001_1"
                 "Sh_Kar_37_2"           "Sh_Lata_078_1"         "Sh_PEM_0103_1"
                 "Sh_PEM_0104_1"         "Sh_PEM_0130_1"         "Sh_Tiag_272_1"
                 "Sh_UNG_0038_1"         "Sh_UNG_0099_1"         "Sh_UNG_0121_1" 
                 "ERR103048"                                                   )

for INDIVIDUAL in "${INDIVIDUALS[@]}"; do

    #use awk to turn into threshold (min 5?)

    #mask low coverage regions of the (targeted) portions of the genome
    
    #extract the cds for each gene

    #combine into a single 
    extract each cds for each gene
    combine cds for each gene










#create one file of phased snps (for all autosomes)
vcfcombine \
    SM_V7_1_beagle_all_samples.vcf \
    SM_V7_2_beagle_all_samples.vcf \
    SM_V7_3_beagle_all_samples.vcf \
    SM_V7_4_beagle_all_samples.vcf \
    SM_V7_5_beagle_all_samples.vcf \
    SM_V7_6_beagle_all_samples.vcf \
    SM_V7_7_beagle_all_samples.vcf \
    >auto_beagle_all_samples.vcf
    #370,770

#now get the phased invadolysin haplotypes
 vcftools \
    --vcf auto_beagle_all_samples.vcf \
    --bed Smp_127030_cds.bed \
    --recode \
    --recode-INFO-all \
    --stdout \
    >Smp_127030_phased.vcf



#now get this vcf for phased invadolysin haplotypes
python ../../scripts/diploid_to_haploid_vcf.py \
    Smp_127030_phased.vcf \
    Smp_127030_phased_hap_A.vcf \
    Smp_127030_phased_hap_B.vcf

#convert each haplotype to a schMan invadolysin vcf file
for HAPLOTYPE in A B; do
    
    #in out files
    IN_VCF="Smp_127030_phased_hap_"$HAPLOTYPE".vcf"
    OUT_VCF="Smp_127030_hap_"$HAPLOTYPE"_schHae_coords.vcf"

    #lift coordinates
    python ../../scripts/schMan_to_schHae_vcf_coords.py $IN_VCF tmp.vcf 

    #strip header
    grep -v "contig=<ID=" tmp.vcf >headerless.vcf

    #add contigs for schHae to header
    ${ENVIRONMENTS["TITAN SINGULARITY"]} \
        gatk SelectVariants \
            -R $HAE_GENOME \
            -V headerless.vcf \
            -O header.vcf

    #sort to generate final vcf file
    ${ENVIRONMENTS["TITAN SINGULARITY"]} \
        gatk SortVcf \
            -R $HAE_GENOME \
            -I header.vcf \
            -O $OUT_VCF

    #clean up
    rm tmp.vcf 
    rm headerless.vcf 
    rm header.vcf*
    rm out.log

done

if [ -f Smp_127030_haplotypes_cds.fasta ]; then
    rm Smp_127030_haplotypes_cds.fasta
fi

for SAMPLE in $(cat ../samples.list); do
    for HAPLOTYPE in A B; do

        #in out files
        IN_VCF="Smp_127030_hap_"$HAPLOTYPE"_schHae_coords.vcf"
        
        #get the vcf
        vcftools \
            --vcf $IN_VCF \
            --indv $SAMPLE \
            --recode \
            --recode-INFO-all \
            --stdout \
            >sample.vcf

        #get the fasta
        gatk \
           -T FastaAlternateReferenceMaker \
           -R KL250964.1_masked.fasta \
           -IUPAC $SAMPLE\
           -o sample.fasta \
           -V sample.vcf

        #clean up fasta for extraction
        sed -i "1c>KL250964.1" sample.fasta
        samtools faidx sample.fasta

        #extract each cds and them combine into one
        bedtools getfasta \
            -bed Smp_127030_cds_schHae.bed \
            -fi sample.fasta \
            -fo exons.fasta

        #combine into a single sequence
         grep -v "^>" exons.fasta \
            | awk 'BEGIN { ORS=""; print" >Sequence_name\n" } { print }' > cds.fasta

        #change the header
        sed -i "1c>$HAPLOTYPE.$SAMPLE" cds.fasta
        
        #add it to combined fasta file
        cat cds.fasta >>Smp_127030_haplotypes_cds.fasta
        echo >>Smp_127030_haplotypes_cds.fasta

        #clean up
        rm sample.fasta
        rm sample.vcf*
        rm exons.fasta*
        rm cds.fasta
        rm out.log
    done
done

#add schMan to the file
cat schMan_Smp_127030_combined_cds.fasta >>Smp_127030_haplotypes_cds.fasta
echo >>Smp_127030_haplotypes_cds.fasta

#get unique haplotypes
fasta_formatter -t -i Smp_127030_haplotypes_cds.fasta -o Smp_127030_haplotypes_cds.tab

sort -u -k2 Smp_127030_haplotypes_cds.tab | awk '{print" >"$1"\n"$2}' >Smp_127030_uniq_haplotypes_cds.fasta

#examined tree and the frequency of sequences to find a reduced subset for PAML 
# runs

grep -e A.ERR310940 \
    -e A.Sh.TZ_PEM0104.1 \
    -e A.ERR103051 \
    -e A.ERR103048 \
    -e A.Sh.NE_Tiag-272.1 \
    -e schMan \
    Smp_127030_haplotypes_cds.tab \
    | awk '{print" >"$1"\n"$2}' \
        >reduced_Smp_127030_uniq_haplotypes_cds.fasta

#find the best tree
raxmlHPC \
    -m GTRCAT \
    -p 12345 \
    -s reduced_Smp_127030_uniq_haplotypes_cds.fasta \
    -# 1000 \
    -n reduced_Smp_127030_uniq_haplotypes_cds \
    -o schMan

#clean up replicates
mkdir raxml_reduced_trees
mv RAxML_*.reduced_Smp_127030_uniq_haplotypes_cds.RUN.* raxml_reduced_trees/


#now start paml/codeml analyses by making selection and control files
# first get the best raxml tree and and indicate branch of interest
cat RAxML_bestTree.reduced_Smp_127030_uniq_haplotypes_cds
#((((A.ERR103051:0.00618220804948685539,A.Sh.TZ_PEM0104.1:0.00280672070595951260):0.00246171314196946042,A.ERR310940:0.00079230556832455656):0.00367147403398323602,(A.Sh.NE_Tiag-272.1:0.00060498047741618143,A.ERR103048:0.00057900965475889772):0.00217271227828123772):0.02950785716451633997,schMan:0.02950785716451633997);


#added branch label for PAML
echo "((((A.ERR103051:0.00618220804948685539,A.Sh.TZ_PEM0104.1:0.00280672070595951260):0.00246171314196946042,A.ERR310940:0.00079230556832455656):0.00367147403398323602,(A.Sh.NE_Tiag-272.1 #1:0.00060498047741618143,A.ERR103048:0.00057900965475889772)#1:0.00217271227828123772):0.02950785716451633997,schMan:0.02950785716451633997);" \
    >reduced_Smp_127030_uniq_haplotypes_cds.tree

# now convert fasta to phylip
fasta_formatter \
    -t \
    -i reduced_Smp_127030_uniq_haplotypes_cds.fasta \
    -o reduced_Smp_127030_uniq_haplotypes_cds.phy

#then add " 6 2262\n" to the beginning to make it a phylip formatted alignment

#convert all Ns to ?s in the sequence
sed -i 's/N/?/g' reduced_Smp_127030_uniq_haplotypes_cds.phy
#then fix a the ?E on the NE sample
sed -i 's/N?/NE/' reduced_Smp_127030_uniq_haplotypes_cds.phy
#fix tabs to spaces
sed -i 's/\t/    /' reduced_Smp_127030_uniq_haplotypes_cds.phy

#make the neutral and selection codeml control files and upload
codeml neutral.ctl
codeml selection.ctl


grep lnL reduced_Smp_127030_uh-branch_site_neutral.H0.mlc 
#lnL(ntime: 10  np: 14):  -2845.725768      +0.000000
grep lnL reduced_Smp_127030_uh-branch_site_selection.H1.mlc
#lnL(ntime: 10  np: 15):  -2845.103222      +0.000000

#LRT = 2*(-2845.103222 - -2845.725768) = 1.245092
chi2
#Chi-square critical values

#                                Significance level

# DF    0.9950   0.9750   0.9000   0.5000   0.1000   0.0500   0.0100   0.0010
#  1    0.0000   0.0010   0.0158   0.4549   2.7055   3.8415   6.6349  10.8276

#reject model with selection 
#now run with ALL uniq haplotypes-----------------------------------------------
#make the tree
raxmlHPC-PTHREADS \
        -T 2 \
        -m GTRCAT \
        -p 12345 \
        -# 1000 \
        -s Smp_127030_uniq_haplotypes_cds.fasta \
        -n Smp_127030_uniq_haplotypes_cds \
        -o schMan

#clean up a bit
mkdir raxml_all_uniq_haplotypes_trees
mv RAxML_*Smp_127030_uniq_haplotypes_cds* raxml_all_uniq_haplotypes_trees/

cat raxml_all_uniq_haplotypes_trees/RAxML_bestTree.Smp_127030_uniq_haplotypes_cds
#((A.ERR310940:0.00133451934108898816,((((A.ERR310937:0.00063928626593180151,A.Sh.NE_NG-011.1:0.00000100000050002909):0.00257661528568875432,A.ERR119612:0.00000100000050002909):0.00064065471620128247,((B.Sh.NE_Dai-013.3:0.00128854550202245474,(A.ERR103048:0.00064235093003460193,((A.Sh.NE_Kar-001.1:0.00064390292987796566,(A.ERR037800:0.00000100000050002909,A.Sh.NE_YK-029.2:0.00064491744825488509):0.00000100000050002909):0.00064390295564639786,A.Sh.NE_Lata-078.1:0.00000100000050002909):0.00000100000050002909):0.00064503434630353213):0.00064247399058681251,B.ERR119612:0.00000100000050002909):0.00000100000050002909):0.00390783835836038505,(A.Sh.TZ_PEM0103.1:0.00063638485263258643,(B.Sh.TZ_PEM0120.1:0.00063722371357710162,((((A.Sh.TZ_UNG0117.1:0.00000100000050002909,(A.Sh.TZ_UNG0139.1:0.00127094925457619102,(B.Sh.TZ_PEM0099.2:0.00000100000050002909,(A.Sh.TZ_PEM0133.1:0.00126873285883879398,((B.ERR103051:0.00126866320274240702,(A.ERR103051:0.00126896197669645774,A.ERR539855:0.00000100000050002909):0.00190872044923837223):0.00000100000050002909,(B.Sh.TZ_UNG0087.2:0.00000100000050002909,(B.Sh.TZ_PEM0103.1:0.00000100000050002909,(B.Sh.TZ_PEM0145.3:0.00064087518460384724,(((((A.Sh.TZ_PEM0099.2:0.00063772546284204606,A.Sh.TZ_UNG0134.1:0.00000100000050002909):0.00063606639400788645,A.Sh.TZ_PEM0139.2:0.00127830321049212498):0.00000100000050002909,(((B.Sh.TZ_UNG0129.2:0.00321759943031953719,(((B.Sh.TZ_PEM0171.1:0.00127245111283949314,A.Sh.TZ_PEM0063.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0111.1:0.00063465392391950383):0.00127404735379609566,B.Sh.TZ_PEM0126.1:0.00000100000050002909):0.00063514446530332769):0.00000100000050002909,A.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0133.1:0.00063519043260941411):0.00063603530497547388):0.00000100000050002909,A.Sh.TZ_PEM0154.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0125.3:0.00192060815928473847):0.00063695393764512983):0.00063467463078119637):0.00063447937791556610):0.00063365040618427440):0.00000100000050002909):0.00063540095162459905):0.00000100000050002909):0.00191099797591265039):0.00449467180776074327,(B.Sh.TZ_UNG0038.1:0.00063499030110272434,((A.Sh.TZ_PEM0079.1:0.00127841087972491228,(B.Sh.TZ_PEM0110.1:0.00516891239923836752,B.Sh.TZ_UNG0111.1:0.00000100000050002909):0.00063494227797514409):0.00000100000050002909,(B.Sh.TZ_PEM0079.1:0.00000100000050002909,(((((B.Sh.TZ_UNG0076.1:0.00063587060918523499,A.Sh.TZ_PEM0094.2:0.00000100000050002909):0.00063717516739843249,(B.Sh.TZ_UNG0117.1:0.00063469068797979773,A.Sh.TZ_PEM0076.1:0.00127420199109140836):0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0139.2:0.00063552392336664316):0.00063733525493169072,A.Sh.TZ_UNG0087.2:0.00063578704845731245):0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00063577085734859925):0.00000100000050002909,B.Sh.TZ_UNG0146.1:0.00000100000050002909):0.00063663384739386401,A.Sh.TZ_PEM0145.3:0.00000100000050002909):0.00063906355174048407):0.00000100000050002909):0.00258673892333877772):0.00125809983233602274):0.03392093211114020901,schMan:0.03392093211114020901);

#add the foreground branch label
echo "((A.ERR310940:0.00133451934108898816,((((A.ERR310937:0.00063928626593180151,A.Sh.NE_NG-011.1:0.00000100000050002909):0.00257661528568875432,A.ERR119612:0.00000100000050002909):0.00064065471620128247,((B.Sh.NE_Dai-013.3:0.00128854550202245474,(A.ERR103048:0.00064235093003460193,((A.Sh.NE_Kar-001.1:0.00064390292987796566,(A.ERR037800:0.00000100000050002909,A.Sh.NE_YK-029.2:0.00064491744825488509):0.00000100000050002909):0.00064390295564639786,A.Sh.NE_Lata-078.1:0.00000100000050002909):0.00000100000050002909):0.00064503434630353213):0.00064247399058681251,B.ERR119612:0.00000100000050002909):0.00000100000050002909)#1:0.00390783835836038505,(A.Sh.TZ_PEM0103.1:0.00063638485263258643,(B.Sh.TZ_PEM0120.1:0.00063722371357710162,((((A.Sh.TZ_UNG0117.1:0.00000100000050002909,(A.Sh.TZ_UNG0139.1:0.00127094925457619102,(B.Sh.TZ_PEM0099.2:0.00000100000050002909,(A.Sh.TZ_PEM0133.1:0.00126873285883879398,((B.ERR103051:0.00126866320274240702,(A.ERR103051:0.00126896197669645774,A.ERR539855:0.00000100000050002909):0.00190872044923837223):0.00000100000050002909,(B.Sh.TZ_UNG0087.2:0.00000100000050002909,(B.Sh.TZ_PEM0103.1:0.00000100000050002909,(B.Sh.TZ_PEM0145.3:0.00064087518460384724,(((((A.Sh.TZ_PEM0099.2:0.00063772546284204606,A.Sh.TZ_UNG0134.1:0.00000100000050002909):0.00063606639400788645,A.Sh.TZ_PEM0139.2:0.00127830321049212498):0.00000100000050002909,(((B.Sh.TZ_UNG0129.2:0.00321759943031953719,(((B.Sh.TZ_PEM0171.1:0.00127245111283949314,A.Sh.TZ_PEM0063.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0111.1:0.00063465392391950383):0.00127404735379609566,B.Sh.TZ_PEM0126.1:0.00000100000050002909):0.00063514446530332769):0.00000100000050002909,A.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0133.1:0.00063519043260941411):0.00063603530497547388):0.00000100000050002909,A.Sh.TZ_PEM0154.1:0.00000100000050002909):0.00000100000050002909,A.Sh.TZ_UNG0125.3:0.00192060815928473847):0.00063695393764512983):0.00063467463078119637):0.00063447937791556610):0.00063365040618427440):0.00000100000050002909):0.00063540095162459905):0.00000100000050002909):0.00191099797591265039):0.00449467180776074327,(B.Sh.TZ_UNG0038.1:0.00063499030110272434,((A.Sh.TZ_PEM0079.1:0.00127841087972491228,(B.Sh.TZ_PEM0110.1:0.00516891239923836752,B.Sh.TZ_UNG0111.1:0.00000100000050002909):0.00063494227797514409):0.00000100000050002909,(B.Sh.TZ_PEM0079.1:0.00000100000050002909,(((((B.Sh.TZ_UNG0076.1:0.00063587060918523499,A.Sh.TZ_PEM0094.2:0.00000100000050002909):0.00063717516739843249,(B.Sh.TZ_UNG0117.1:0.00063469068797979773,A.Sh.TZ_PEM0076.1:0.00127420199109140836):0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0114.3:0.00000100000050002909):0.00000100000050002909,B.Sh.TZ_PEM0139.2:0.00063552392336664316):0.00063733525493169072,A.Sh.TZ_UNG0087.2:0.00063578704845731245):0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00063577085734859925):0.00000100000050002909,B.Sh.TZ_UNG0146.1:0.00000100000050002909):0.00063663384739386401,A.Sh.TZ_PEM0145.3:0.00000100000050002909):0.00063906355174048407):0.00000100000050002909):0.00258673892333877772):0.00125809983233602274):0.03392093211114020901,schMan:0.03392093211114020901);" \
    >Smp_127030_uniq_haplotypes_cds.tree

#make the phylip file
fasta_formatter \
    -t \
    -i Smp_127030_uniq_haplotypes_cds.fasta \
    -o Smp_127030_uniq_haplotypes_cds.tab

cut -f1 Smp_127030_uniq_haplotypes_cds.tab >names
cut -f2 Smp_127030_uniq_haplotypes_cds.tab | sed 's/N/\?/g' >seqs

paste names seqs >Smp_127030_uniq_haplotypes_cds.phy

#then add " 50 2262\n" to the beginning to make it a phylip formatted alignment

#fix tabs to spaces
sed -i 's/\t/    /' Smp_127030_uniq_haplotypes_cds.phy

#make the ctl files
sed 's/reduced_//' \
    <neutral.ctl \
    >neutral_all_uniq_haps/neutral_all_uniq_haps.ctl
cp Smp_127030_uniq_haplotypes_cds.tree neutral_all_uniq_haps/
cp Smp_127030_uniq_haplotypes_cds.phy neutral_all_uniq_haps/

sed 's/reduced_//' \
    <selection.ctl \
    >selection_all_uniq_haps/selection_all_uniq_haps.ctl
cp Smp_127030_uniq_haplotypes_cds.tree selection_all_uniq_haps/
cp Smp_127030_uniq_haplotypes_cds.phy selection_all_uniq_haps/

sed 's/reduced_//' \
    <dnds.ctl \
    >dnds_all_uniq_haps/dnds_all_uniq_haps.ctl    
cp Smp_127030_uniq_haplotypes_cds.tree dnds_all_uniq_haps/
cp Smp_127030_uniq_haplotypes_cds.phy dnds_all_uniq_haps/

codeml neutral_all_uniq_haps.ctl
codeml selection_all_uniq_haps.ctl
codeml dnds_all_uniq_haps.ctl

#site-test
#neutral lnL(ntime: 98  np:102):  -3392.298001      +0.000000
#selection lnL(ntime: 98  np:103):  -3392.298006      +0.000000

#SM_V7_4 20033013        KL250964.1:40108 
