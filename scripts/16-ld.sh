#clean and process reads to the haematobium genome
source /master/nplatt/schisto_hybridization/scripts/set_env.sh
source activate snp_calling

cd $RESULTS_DIR

mkdir ld

cd ld

vcftools \
    --vcf ../build_snp_panel/auto_maf.vcf \
    --maf 0.05 \
    --recode \
    --stdout \
    >auto_maf05.vcf

vcftools \
    --vcf SM_V7_1_maf05.vcf \
    --keep tz.samples \
    --maf 0.1 \
    --recode \
    --stdout \
    >tz_SM_V7_1_maf10.vcf


plink \
    --vcf tz_SM_V7_1_maf10.vcf \
    --out tz_SM_V7_1_maf10 \
    --recode12 \
    --allow-extra-chr

plink \
    --r2 \
    --file tz_SM_V7_1_maf10 \
    --out tz_SM_V7_1_maf10 \
    --allow-extra-chr \
    --ld-window-r2 0 \
    --ld-window 22000 \
    --ld-window-kb 88882

R
library(scales)
test<-read.table("SM_V7_1_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_1_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,8.9e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_2_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_2_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,4.9e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_3_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_3_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,5.1e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_4_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_4_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,4.8e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_5_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_5_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,2.6e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_6_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_6_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,2.5e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()


test<-read.table("SM_V7_7_tz_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_7_tz_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,19e6), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

#ne
test<-read.table("SM_V7_1_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_1_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,8.9e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_2_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_2_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,4.9e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_3_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_3_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,5.1e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_4_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_4_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,4.8e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_5_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_5_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,2.6e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_6_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_6_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,2.5e7), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()

test<-read.table("SM_V7_7_ne_haem.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('SM_V7_7_ne_haem.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,19e6), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
dev.off()


vcftools \
    --vcf auto_maf05.vcf \
    --bed start.bed \
    --recode-INFO GT \
    --recode \
    --stdout \
    >start.vcf

vcftools \
    --vcf auto_maf05.vcf \
    --maf 0.1 \
    --bed end.bed \
    --recode \
    --stdout \
    
    >end.vcf


#fitting
awk '{ if ($7!=0) print $0}' SM_V7_7_tz_haem.ld >notZero.tbl

R
tz_7<-read.table("notZero.tbl", header=TRUE)
dist<-abs(tz_7$BP_A - tz_7$BP_B) 
exponential.model<-lm(log(tz_7$R2) ~ dist)
#need to get rid of R2=0

dist<- seq(1, 19e6, 0.5e6)
counts.exponential2 <- exp(predict(exponential.model,list(dist=distvalues)))
png("test.png")
plot(x=dist, y=tz_7$R2, xlim=c(0,19e6), ylim=c(0,1), cex=0.1, pch=19, col=alpha("black", 0.1))
lines(distvalues, counts.exponential2,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")
dev.off()

plot(hexbin(x=test_dist, y=test$R2, xbins = 50))



mkdir schMan_ld
for i in $(seq 1 7); do
    for POP in tz ne; do
    CHR=SM_V7_"$i"

        #get vcf chromosome
        vcftools \
            --vcf ../build_snp_panel/auto_maf.vcf \
            --keep $POP.samples \
            --maf 0.1 \
            --chr $CHR \
            --recode \
            --stdout \
            >$CHR"_"$POP".vcf"

        #convert to ped
        plink \
            --threads 6 \
            --vcf $CHR"_"$POP".vcf" \
            --out $CHR"_"$POP \
            --recode12 \
            --allow-extra-chr

        #calc R2 for all snps on the chr
        plink \
            --threads 6 \
            --r2 \
            --file $CHR"_"$POP \
            --out $CHR"_"$POP \
            --allow-extra-chr \
            --ld-window-r2 0 \
            --ld-window 75000 \
            --ld-window-kb 100000


        awk '{print $0"\t"$5-$2}' $CHR"_"$POP".ld" >schMan_ld/$CHR"_"$POP".ld"
    done
done



#convert vcf back to SchHaem coordinates
grep -v "#" ../build_snp_panel/auto_maf.vcf \
    | awk  '{print $3"\t"$0}' \
    | sed 's/:/\t/' \
    | awk  '{$3=$3":"$4;; print $0}' \
    | sed 's/ /\t/g' \
    | cut --complement -f4,5 \
    >vcf.entries

grep "#" ../build_snp_panel/auto_maf.vcf | grep -v "contig=<ID=" >vcf.sm_header   

cat vcf.sm_header vcf.entries >broke_header.vcf

#add contigs for sman to header
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SelectVariants \
        -R $HAE_GENOME \
        -V broke_header.vcf \
        -O header.vcf
        #75,704 loci remaining

#now sort
${ENVIRONMENTS["SINGULARITY"]} \
    gatk SortVcf \
        -R $HAE_GENOME \
        -I header.vcf \
        -O auto_maf_schHae.vcf
        #75,704 snp loci across all (sman homologous) chromosomes

#create list of contigs
grep ">" $HAE_GENOME | sed 's/>//' | cut -f1 -d" " >schHae_contigs.list

mkdir schHae_ld
#calc ld for all schHae contigs
for CHR in $(cat schHae_contigs.list); do

    #get vcf chromosome
    vcftools \
        --vcf auto_maf_schHae.vcf \
        --keep all.samples \
        --chr $CHR \
        --recode \
        --stdout \
        >$CHR".vcf"

    #convert to ped
    plink \
        --threads 6 \
        --vcf $CHR".vcf" \
        --out $CHR \
        --recode12 \
        --allow-extra-chr

    #calc R2 for all snps on the chr
    plink \
        --threads 6 \
        --r2 \
        --file $CHR \
        --out $CHR \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000

    #calc distance
    awk '{print $0"\t"$5-$2}' $CHR.ld >schHae_ld/$CHR.ld
    rm $CHR.ped $CHR.map $CHR.nosex $CHR.log $CHR.vcf $CHR.ld

done

rm *temp*

#when finished
cat schHae_ld/*.ld | awk '{print $7","$8}' | grep -v "R2" >schHae_ld.csv
tar -czf schHae_ld.tgz schHae_ld/
rm -r schHae_ld/
#produces a rather small ~75MB file

for i in $(seq 1 7); do

    CHR=SM_V7_"$i"

    awk '{print $7","$8 }' $CHR"_haem.ld" | grep -v R2 >$CHR"_haem.ld.csv" &

done
wait


##if downsampling in python
#from random import random
#out_file = open('niger_snps_schMan_chr1_downsampled', 'w')

#lines = [line for line in open("niger_snps_schMan_chr1.ld") if random() >= .005]

#for item in lines:
#  out_file.write("%s\n" % item)

#download to local computer for analysis in R

setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/ld")

niger_ld<-read.table("niger_snps_schMan_autosomal_panel.ld", header=TRUE)
tz_ld<-read.table("tz_snps_schMan_autosomal_panel.ld", header=TRUE)

tz_chr1_ld<-read.table("tz_chr1_snps_schMan_autosomal_panel.ld", header=TRUE)


niger_dist<-abs(niger_ld$BP_A - niger_ld$BP_B) 

tz_dist<-abs(tz_ld$BP_A - tz_ld$BP_B) 

tz_chr1_dist<-abs(tz_chr1_ld$BP_A - tz_chr1_ld$BP_B) 

max(tz_dist)
max(tz_chr1_dist)

plot(x=tz_chr1_dist, y=tz_chr1_ld$R2, xlim=c(0,1000000), ylim=c(0,1))





setwd("C:/Users/nplatt/Dropbox/work/projects/sH_hybridization/results/ld")

haem_ld<-read.csv("schHae_ld.csv", header=FALSE, sep=",")
png('schHae_ld.png')
plot(hexbin(x=haem_ld$V2, y=haem_ld$V1, xbins = 50))
dev.off()

man_ld<-read.table("schMan_ld.csv", header=FALSE, sep=",")
png('schMan_ld.png')
plot(hexbin(x=man_ld$V2, y=man_ld$V1, xbins = 50))
dev.off()


${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        --select-random-fraction 0.1 \
        -V SM_V7_1_haem.vcf \
        -O SM_V7_1_haem_10percSub.vcf

    vcftools \
        --vcf SM_V7_1_haem_10percSub.vcf \
        --keep ../niger.samples \
        --recode \
        --stdout \
        >niger.vcf

    vcftools \
        --vcf SM_V7_1_haem_10percSub.vcf \
        --keep ../tz.samples \
        --recode \
        --stdout \
        >tz.vcf

plink \
        --threads 6 \
        --vcf tz.vcf \
        --out tz \
        --recode12 \
        --allow-extra-chr

plink \
        --threads 6 \
        --vcf niger.vcf \
        --out niger \
        --recode12 \
        --allow-extra-chr

    #calc R2 for all snps on the chr
    plink \
        --threads 6 \
        --r2 \
        --file tz \
        --out tz \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000

    plink \
        --threads 6 \
        --r2 \
        --file niger \
        --out niger \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000

awk '{print $7","$5-$2 }' niger.ld | grep -v R2 >niger.csv &
awk '{print $7","$5-$2 }' tz.ld | grep -v R2 >tz.csv &


cp ../../beagle/SM_V7_1_beagle.vcf .

${ENVIRONMENTS["TITAN SINGULARITY"]} \
    gatk SelectVariants \
        --select-random-fraction 0.1 \
        -V SM_V7_1_beagle.vcf \
        -O beagle10.vcf

vcftools \
        --vcf beagle10.vcf \
        --keep ../niger.samples \
        --recode \
        --stdout \
        >beagle_ne.vcf

    vcftools \
        --vcf beagle10.vcf \
        --keep ../tz.samples \
        --recode \
        --stdout \
        >beagle_tz.vcf

plink \
        --threads 6 \
        --vcf beagle_tz.vcf \
        --out beagle_tz \
        --recode12 \
        --allow-extra-chr

plink \
        --threads 6 \
        --vcf beagle_ne.vcf \
        --out beagle_ne \
        --recode12 \
        --allow-extra-chr


    plink \
        --threads 6 \
        --r2 \
        --file beagle_ne \
        --out beagle_ne \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000


    plink \
        --threads 6 \
        --r2 \
        --file beagle_tz \
        --out beagle_tz \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000

awk '{print $7","$5-$2 }' beagle_tz.ld | grep -v R2 >beagle_tz.csv &
awk '{print $7","$5-$2 }' beagle_ne.ld | grep -v R2 >beagle_ne.csv &


cat SM_V7_1_haem.ld | awk '{if ($7>0.6 && $5-$2>250000) print $0}'
vcftools --vcf ../beagle/SM_V7_1_beagle.vcf --keep tz.samples --maf 0.05 --recode --stdout >test.vcf

plink \
        --threads 6 \
        --vcf test.vcf \
        --out test \
        --recode12 \
        --allow-extra-chr


    plink \
        --threads 6 \
        --r2 \
        --file test \
        --out test \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000


vcftools --vcf ../beagle/SM_V7_1_beagle.vcf --keep niger.samples --maf 0.05 --recode --stdout >test_ne.vcf

plink \
        --threads 6 \
        --vcf test_ne.vcf \
        --out test_ne \
        --recode12 \
        --allow-extra-chr


    plink \
        --threads 6 \
        --r2 \
        --file test \
        --out test \
        --allow-extra-chr \
        --ld-window-r2 0 \
        --ld-window 75000 \
        --ld-window-kb 100000


test<-read.table("test.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('test_tz_chr1_nobeagle.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()



#code that recreates LD decay with no weird linkage
#vcftools \
#    --vcf ../beagle_impute/SM_V7_1_beagle.vcf \
#    --maf 0.05 \
#    --keep tz.samples \
#    --recode \
#    --stdout \
#    >tz_beagle_SM_V7_1.vcf

#plink \
#    --vcf tz_beagle_SM_V7_1.vcf \
#    --out tz_beagle_SM_V7_1 \
#    --recode12 \
#    --allow-extra-chr

#plink \
#    --r2 \
#    --file tz_beagle_SM_V7_1 \
#    --out tz_beagle_SM_V7_1 \
#    --allow-extra-chr \
#    --ld-window-r2 0 \
#    --ld-window 22000 \
#    --ld-window-kb 88882

#vcftools     --vcf ../../beagle/old/SM_V7_1_beagle.vcf     --maf 0.05     --keep tz.samples     --recode     --stdout     >test.vcf
#plink     --vcf test.vcf     --out test     --recode12     --allow-extra-chr
#plink     --r2     --file test     --out test     --allow-extra-chr     --ld-window-r2 0     --ld-window 22000     --ld-window-kb 88882



#test<-read.table("test.ld", header=TRUE)
#test_dist<-abs(test$BP_A - test$BP_B) 
#png('test_tz_chr1_beagle.ld.png')
#plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
#dev.off()


vcftools \
    --vcf ../beagle/SM_V7_1_beagle.vcf \
    --maf 0.05 \
    --keep tz.samples \
    --recode \
    --stdout \
    >tz_beagle_SM_V7_1_new.vcf

plink \
    --vcf tz_beagle_SM_V7_1_new.vcf \
    --out tz_beagle_SM_V7_1_new \
    --recode12 \
    --allow-extra-chr

plink \
    --r2 \
    --file tz_beagle_SM_V7_1_new \
    --out tz_beagle_SM_V7_1_new \
    --allow-extra-chr \
    --ld-window-r2 0 \
    --ld-window 22000 \
    --ld-window-kb 88882



test<-read.table("test.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('tz_beagle_SM_V7_1_new.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()


vcftools \
    --vcf ../../beagle/SM_V7_1_beagle.vcf \
    --snps old.list \
    --keep tz.samples \
    --recode \
    --stdout \
    >new-vcf_old-list.vcf

plink \
    --vcf new-vcf_old-list.vcf \
    --out new-vcf_old-list \
    --recode12 \
    --allow-extra-chr

plink \
    --r2 \
    --file new-vcf_old-list \
    --out new-vcf_old-list \
    --allow-extra-chr \
    --ld-window-r2 0 \
    --ld-window 22000 \
    --ld-window-kb 88882

test<-read.table("new-vcf_old-list.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('new-vcf_old-list.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()



vcftools \
    --vcf ../../beagle/SM_V7_1_beagle.vcf \
    --maf 0.1 \
    --keep tz.samples \
    --recode \
    --stdout \
    >tz_beagle_SM_V7_1_maf10.vcf

plink \
    --vcf tz_beagle_SM_V7_1_maf10.vcf \
    --out tz_beagle_SM_V7_1_maf10 \
    --recode12 \
    --allow-extra-chr

plink \
    --r2 \
    --file tz_beagle_SM_V7_1_maf10 \
    --out tz_beagle_SM_V7_1_maf10 \
    --allow-extra-chr \
    --ld-window-r2 0 \
    --ld-window 22000 \
    --ld-window-kb 88882



test<-read.table("tz_beagle_SM_V7_1_maf10.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('tz_beagle_SM_V7_1_maf10.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()


vcftools \
    --vcf ../../beagle/SM_V7_1_beagle.vcf \
    --maf 0.05 \
    --recode \
    --stdout \
    >beagle_SM_V7_1_maf05.vcf

vcftools \
    --vcf beagle_SM_V7_1_maf05.vcf \
    --keep tz.samples \
    --maf 0.05 \
    --recode \
    --stdout \
    >tz_beagle_SM_V7_1_maf05.vcf


plink \
    --vcf tz_beagle_SM_V7_1_maf05.vcf \
    --out tz_beagle_SM_V7_1_maf05 \
    --recode12 \
    --allow-extra-chr

plink \
    --r2 \
    --file tz_beagle_SM_V7_1_maf05 \
    --out tz_beagle_SM_V7_1_maf05 \
    --allow-extra-chr \
    --ld-window-r2 0 \
    --ld-window 22000 \
    --ld-window-kb 88882


R
test<-read.table("tz_beagle_SM_V7_1_maf05.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('tz_beagle_SM_V7_1_maf05.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()




vcftools \
    --vcf ../../beagle/SM_V7_1_beagle.vcf \
    --maf 0.05 \
    --recode \
    --stdout \
    >beagle_SM_V7_1_maf05.vcf

vcftools \
    --vcf beagle_SM_V7_1_maf05.vcf \
    --keep niger.samples \
    --maf 0.05 \
    --recode \
    --stdout \
    >ne_beagle_SM_V7_1_maf05.vcf


plink \
    --vcf ne_beagle_SM_V7_1_maf05.vcf \
    --out ne_beagle_SM_V7_1_maf05 \
    --recode12 \
    --allow-extra-chr

plink \
    --r2 \
    --file ne_beagle_SM_V7_1_maf05 \
    --out ne_beagle_SM_V7_1_maf05 \
    --allow-extra-chr \
    --ld-window-r2 0 \
    --ld-window 22000 \
    --ld-window-kb 88882


R
test<-read.table("ne_beagle_SM_V7_1_maf05.ld", header=TRUE)
test_dist<-abs(test$BP_A - test$BP_B) 
png('ne_beagle_SM_V7_1_maf05.ld.png')
plot(x=test_dist, y=test$R2, xlim=c(0,10e7), ylim=c(0,1))
dev.off()



