#read in freq tables from vcftools
#clean up
#get rid of non-abba baba sites
#calculate abba baba scores
#calculate genome wide D
#jacknife
#calculate SD, P, Z values
import numpy
import scipy.stats

txt_in=open("bov_abba_babba", 'r')

chr=[]
pos=[]
tz_freq=[]
ne_freq=[]
bo_freq=[]
ma_freq=[]
abba=[]
baba=[]

for site_entry in txt_in:
       
    site_entry = site_entry.rstrip()
    
    (site_chr, site_pos, tanz, nigr, bovi, marg)=site_entry.split("\t")
        
    site_chr=str(site_chr)
    site_pos=int(site_pos)
    tanz=float(tanz)
    nigr=float(nigr)
    bovi=float(bovi)
    marg=float(marg)
    
    site_abba=(1-tanz)*nigr*bovi
    site_baba=tanz*(1-nigr)*bovi
    
    #append to list
    chr.append(site_chr)
    pos.append(site_pos)
    tz_freq.append(tanz)
    ne_freq.append(nigr)
    bo_freq.append(bovi)
    ma_freq.append(marg)
    abba.append(site_abba)
    baba.append(site_baba)

#calculate genome wide D    
d_stat=(sum(abba)-sum(baba))/(sum(abba)+sum(baba))

#generate files excluding 
chr_length={}
chr_length["SM_V7_1"]=88881357
chr_length["SM_V7_2"]=48130368
chr_length["SM_V7_3"]=50458499
chr_length["SM_V7_4"]=47279781
chr_length["SM_V7_5"]=25256119
chr_length["SM_V7_6"]=24989083
chr_length["SM_V7_7"]=19288021

#build jacknife windows to exclude
window_size=3000000
current_position=0

jacknife_chr=[]
jacknife_window_start=[]
jacknife_window_end=[]

for chrom in "SM_V7_1", "SM_V7_2", "SM_V7_3", "SM_V7_4", "SM_V7_5", "SM_V7_6", "SM_V7_7": 
    window_start=0
    window_end=window_start+window_size
    
    while window_end < chr_length[chrom]:
        #print(chrom, window_start, window_end)
        jacknife_chr.append(chrom)     
        jacknife_window_start.append(window_start)
        jacknife_window_end.append(window_end)
        
        window_start=window_end
        window_end=window_end+window_size
                
        if window_end > chr_length[chrom]:
            window_end=chr_length[chrom]
            
            jacknife_chr.append(chrom)     
            jacknife_window_start.append(window_start)
            jacknife_window_end.append(window_end)


#now i have my list of windows to exclude - calculate D for all except in this window

jacknife_out=open("jacknife.csv", 'w')

window_abba=[]
window_baba=[]
jacknife_d=[]

#for each window
for window in range(len(jacknife_chr)):
    window_chr=jacknife_chr[window]
    window_start=jacknife_window_start[window]
    window_end=jacknife_window_end[window]
    window_n=0
    
    #for each sample that falls outside of the window    
    for sample in range(len(tz_freq)):
        
        beyond_window=not(pos[sample] >= window_start and pos[sample] <= window_end)
        diff_chr=not(chr[sample] == window_chr)
        
        if diff_chr or beyond_window:
            window_abba.append(abba[sample])
            window_baba.append(baba[sample])
            window_n=window_n+1
                
        
    window_d=(sum(window_abba)-sum(window_baba))/(sum(window_abba)+sum(window_baba))
    window_entry=str(window_chr) + "," + str(window_start) + "," +  str(window_end) + "," +  str(window_n) + "," + str(window_d) + "\n"
    jacknife_out.write(window_entry)
    jacknife_d.append(window_d)

jacknife_out.close()

#d stats in file
#calculate stdDev, stdErr, Z, P
D_sd=numpy.std(jacknife_d)
D_err=D_sd/numpy.sqrt(len(jacknife_d))
D_Z=d_stat/D_err
D_p <- 2*pnorm(-abs(D_Z)) #haven't done in python but in R p=0 (significant deviation)




#now do dstats on the across the genome
#create windows of X bp
#calculate D for that window
#slide to next window    
window_size=50000
step=10000

slide_window_chr=[]
slide_window_start=[]
slide_window_end=[]

      
for chrom in "SM_V7_1", "SM_V7_2", "SM_V7_3", "SM_V7_4", "SM_V7_5", "SM_V7_6", "SM_V7_7": 
    window_start=0
    window_end=window_start+window_size
    
    while window_end < chr_length[chrom]:
        #print(chrom, window_start, window_end)
        slide_window_chr.append(chrom)     
        slide_window_start.append(window_start)
        slide_window_end.append(window_end)
        
        window_start=window_start+step
        window_end=window_start+window_size
                
        if window_end > chr_length[chrom]:
            window_end=chr_length[chrom]
            
            slide_window_chr.append(chrom)     
            slide_window_start.append(window_start)
            slide_window_end.append(window_end)


#now calculate d for each window
slide_d_out=open("slide_d.csv", 'w')

window_abba=[]
window_baba=[]
slide_d=[]

#for each window
for window in range(len(slide_window_chr)):
    window_chr=slide_window_chr[window]
    window_start=slide_window_start[window]
    window_end=slide_window_end[window]
    window_n=0
    
    #for each sample that falls outside of the window    
    for sample in range(len(tz_freq)):
        
        in_window=pos[sample] >= window_start and pos[sample] <= window_end
        same_chr=chr[sample] == window_chr
        
        if same_chr and in_window:
            window_abba.append(abba[sample])
            window_baba.append(baba[sample])
            window_n=window_n+1
            
    if window_n > 0: 
        window_d=(sum(window_abba)-sum(window_baba))/(sum(window_abba)+sum(window_baba))
        window_entry=str(window_chr) + "," + str(window_start) + "," +  str(window_end) + "," +  str(window_n) + "," + str(sum(window_abba)) + "," + str(sum(window_baba)) + "," + str(window_d)+ "\n"
        
    else:
        window_entry=str(window_chr) + "," + str(window_start) + "," +  str(window_end) + "," +  str(window_n) + "," + "0" + "," + "0" + "\n"
    
    slide_d_out.write(window_entry)
    slide_d.append(window_d)


slide_d_out.close()











