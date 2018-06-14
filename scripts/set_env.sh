#TODO: update all code onto github


#scripts used:
#1 - project setup - makes major dirs and access private data
#2 - get_sra_data - gets data from sra
#3 - filter_and_map_reads_to_haem.sh - maps reads from all species to haematobium
#       and runs reads through trimmomatic
#4 - BQSR-r1 - base quality recal round 1
#5 - BQSR-r2 - base quality recal round 2
#6 - HQ_genotype - gets HQ genotypes from SH samples
#7 - genotype_all_samples - genotypes SH and all outgroup samples at SNP panel
#       and polishes

#~~~~~~~~~~data exploration
#a - admixture
#b - pca
#c - fst and het
#d - phylogeny
#e - 

#set major directories
WORK_DIR=/master/nplatt/schisto_hybridization
RESULTS_DIR=$WORK_DIR/results
LOGS_DIR=$RESULTS_DIR/logs
SCRIPTS_DIR=$RESULTS_DIR/scripts
FILTER_DIR=$RESULTS_DIR/processed_reads/filtered_reads
MAP_DIR=$RESULTS_DIR/processed_reads/sh_mapped_reads
SNP_DIR=$RESULTS_DIR/call_snps
SNP_PANEL_DIR=$RESULTS_DIR/snp_panel


DATA_DIR=$WORK_DIR/data
SEQ_DIR=$DATA_DIR/sequence_data
GENOME_DIR=$DATA_DIR/genome

#set project specific/sample parameters
#set environment and submission variables
declare -A ENVIRONMENTS
declare -A RG_CELL
declare -A RG_LANE
declare -A RG_INDEX
declare -A GENOMES

#base qsub command
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y"    

#various conda/singularity containers for analyses
ENVIRONMENTS=(  ["SINGULARITY"]="singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img"
                ["CONDA_SNPS"]="source activate snp_calling;"
                ["CONDA_HYBRID_ID"]="source activate hybrid_id;"
                ["TITAN SINGULARITY"]="/opt/projects/singularity-2.4.2/bin/singularity exec $WORK_DIR/config/snpCalling_v0.0.8.img"
             )

#species IDs
#ERR119622 bovis -LR - ERS094852 - 
#ERR103048 bovis - ERS074976 - 15.4
#ERR539853 bovis - ERS427204 - 7.9

#ERR310937 curassoni
#ERR119623 curassoni

#ERR084970 haematobium
#ERR037800 haematobium
#SRR433865 haematobium

#ERR539854 intercalatum 
#ERR539856 intercalatum
#ERR119613 intercalatum 

#ERR119612 guineensis
#ERR539850 guineensis
#ERR539852 guineensis

#ERR539851 mattheei
#ERR539855 mattheei
#ERR539857 mattheei
#ERR103051 mattheei - ERS074983  - 20.1

#ERR310940 margrebowiei

#genomes for each species
BOV_GENOME=$GENOME_DIR"/schBov_v1.fa"
CUR_GENOME=$GENOME_DIR"/schCur_v1.fa"
MAN_GENOME=$GENOME_DIR"/schMan_v7.fa"
HAE_GENOME=$GENOME_DIR"/schHae_v1.fa"
GUI_GENOME=$GENOME_DIR"/schGui_v1.fa"
INT_GENOME=$GENOME_DIR"/schInt_v1.fa"
MAR_GENOME=$GENOME_DIR"/schMar_v1.fa"
MAT_GENOME=$GENOME_DIR"/schMat_v1.fa"

#MAN_LINK=$GENOME_DIR"/schMan_v7.fa"
#HAE_LINK=$GENOME_DIR"/schHae_v1.fa"
GENOME_LINKS=(
    ["BOV"]="ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_bovis_v1.0.fa.gz"
    ["CUR"]="ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_curassoni_v1.0_4.fa.gz"
    ["GUI"]="ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_guineensis_v1.0.fa.gz"
    ["INT"]="ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_intercalatum_v1.0.fa.gz"
    ["MAR"]="ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_margrebowiei_v1.5_4.fa.gz"
    ["MAT"]="ftp://ftp.sanger.ac.uk/pub/project/pathogens/HGI/Schistosoma_mattheei_v1.0_4.fa.gz"
    )

SAMPLES=(   "Sm.BR_0447.1"          "Sm.BR_1278.1"          "Sm.BR_2039.1"      
            "ERR103048"             "ERR119622"             "ERR119623" 
            "ERR310937"             "SRR433865"             "ERR037800"
            "ERR084970"             "Sh.NE_Dai-002.1"       "Sh.NE_Dai-010.1" 
            "Sh.NE_Dai-013.3"       "Sh.NE_Dai-031.1"       "Sh.NE_Dai-033.1"
            "Sh.NE_Dai-044.1"       "Sh.NE_Dai-045.1"       "Sh.NE_Dai-051.1"   
            "Sh.NE_Dai-074.1"       "Sh.NE_Dai-146.1"       "Sh.NE_DaiCP-233.1"
            "Sh.NE_DaiCP-276.1"     "Sh.NE_DaiCP-277.2"     "Sh.NE_Doki-029.1" 
            "Sh.NE_Kar-001.1"       "Sh.NE_Kar-002.1"       "Sh.NE_Kar-076.1" 
            "Sh.NE_Kar-096.2"       "Sh.NE_Kar-241.1"       "Sh.NE_Kar-241.2" 
            "Sh.NE_Kar-281.1"       "Sh.NE_Kar-37.2"        "Sh.NE_Lata-007.3" 
            "Sh.NE_Lata-033.1"      "Sh.NE_Lata-078.1"      "Sh.NE_Lata-253.1" 
            "Sh.NE_Lata-275.2"      "Sh.NE_Lata-293.1"      "Sh.NE_Lata-294.1" 
            "Sh.NE_LibTB-009.2"     "Sh.NE_LibTB-010.1"     "Sh.NE_LibTB-022.1" 
            "Sh.NE_LibTB-028.1"     "Sh.NE_LibTB-031.1"     "Sh.NE_NG-011.1" 
            "Sh.NE_NG-06.2"         "Sh.NE_NG-089.1"        "Sh.NE_NG-236.1" 
            "Sh.NE_Seb-076.1"       "Sh.NE_Seb-078.2"       "Sh.NE_Seb-081.2" 
            "Sh.NE_Tiag-272.1"      "Sh.NE_YK-029.2"        "Sh.NE_YK-069.1" 
            "Sh.NE_YK-099.2"        "Sh.NE_YK-248.2"        "Sh.NE_Youri-069.2" 
            "Sh.NE_Youri-091.3"     "Sh.TZ_PEM0063.1"       "Sh.TZ_PEM0075.1" 
            "Sh.TZ_PEM0076.1"       "Sh.TZ_PEM0079.1"       "Sh.TZ_PEM0089.2" 
            "Sh.TZ_PEM0094.2"       "Sh.TZ_PEM0099.2"       "Sh.TZ_PEM0103.1" 
            "Sh.TZ_PEM0104.1"       "Sh.TZ_PEM0106.2"       "Sh.TZ_PEM0108.1" 
            "Sh.TZ_PEM0110.1"       "Sh.TZ_PEM0114.3"       "Sh.TZ_PEM0115.4" 
            "Sh.TZ_PEM0120.1"       "Sh.TZ_PEM0125.1"       "Sh.TZ_PEM0126.1" 
            "Sh.TZ_PEM0127.1"       "Sh.TZ_PEM0128.1"       "Sh.TZ_PEM0130.1" 
            "Sh.TZ_PEM0133.1"       "Sh.TZ_PEM0139.2"       "Sh.TZ_PEM0145.3" 
            "Sh.TZ_PEM0154.1"       "Sh.TZ_PEM0157.3"       "Sh.TZ_PEM0166.1" 
            "Sh.TZ_PEM0171.1"       "Sh.TZ_UNG0006.1"       "Sh.TZ_UNG0038.1"  
            "Sh.TZ_UNG0076.1"       "Sh.TZ_UNG0077.1"       "Sh.TZ_UNG0078.1" 
            "Sh.TZ_UNG0087.2"       "Sh.TZ_UNG0089.3"       "Sh.TZ_UNG0092.3" 
            "Sh.TZ_UNG0099.1"       "Sh.TZ_UNG0102.1"       "Sh.TZ_UNG0111.1" 
            "Sh.TZ_UNG0117.1"       "Sh.TZ_UNG0121.1"       "Sh.TZ_UNG0125.3" 
            "Sh.TZ_UNG0127.1"       "Sh.TZ_UNG0129.2"       "Sh.TZ_UNG0134.1" 
            "Sh.TZ_UNG0137.3"       "Sh.TZ_UNG0139.1"       "Sh.TZ_UNG0142.2" 
            "Sh.TZ_UNG0146.1"       "ERR103051"             "ERR119612" 
            "ERR119613"             "ERR310940"             "ERR539850" 
            "ERR539851"             "ERR539852"             "ERR539853" 
            "ERR539854"             "ERR539855"             "ERR539856" 
            "ERR539857"
        )

SRA_SAMPLES=(   "ERR119622" "ERR103048" "ERR119622" "ERR310937" 
                "ERR119623" "ERR084970" "ERR037800" "SRR433865" 
                "ERR103051" "ERR119612" "ERR119613" "ERR310940" 
                "ERR539850" "ERR539851" "ERR539852" "ERR539853" 
                "ERR539854" "ERR539855" "ERR539856" "ERR539857"
            )

RG_CELL=(  ["ERR119622"]="4"                   ["ERR103048"]="3" 
           ["ERR119623"]="5"                   ["Sm.BR_1278.1"]="C3E5HACXX" 
           ["Sm.BR_0447.1"]="C3E5HACXX"        ["Sm.BR_2039.1"]="C3E5HACXX" 
           ["ERR310937"]="HS27_7296"           ["ERR037800"]="1"
           ["ERR084970"]="2"                   ["SRR433865"]="6"
           ["Sh.NE_Dai-002.1"]="C78YTACXX"     ["Sh.NE_Dai-010.1"]="C78YTACXX"
           ["Sh.NE_Dai-013.3"]="C70LCACXX"     ["Sh.NE_Dai-031.1"]="C78YTACXX"
           ["Sh.NE_Dai-033.1"]="C78YTACXX"     ["Sh.NE_Dai-044.1"]="C70LCACXX"
           ["Sh.NE_Dai-045.1"]="C78YTACXX"     ["Sh.NE_Dai-051.1"]="C78YTACXX"
           ["Sh.NE_Dai-074.1"]="C78YTACXX"     ["Sh.NE_Dai-146.1"]="C78YTACXX"
           ["Sh.NE_DaiCP-233.1"]="C70LCACXX"   ["Sh.NE_DaiCP-276.1"]="C78YTACXX"
           ["Sh.NE_DaiCP-277.2"]="C78YTACXX"   ["Sh.NE_Doki-029.1"]="C78YTACXX"
           ["Sh.NE_Kar-001.1"]="C78YTACXX"     ["Sh.NE_Kar-002.1"]="C78YTACXX"
           ["Sh.NE_Kar-076.1"]="C78YTACXX"     ["Sh.NE_Kar-096.2"]="C78YTACXX"
           ["Sh.NE_Kar-241.1"]="C78YTACXX"     ["Sh.NE_Kar-241.2"]="C78YTACXX"
           ["Sh.NE_Kar-281.1"]="C78YTACXX"     ["Sh.NE_Kar-37.2"]="C70LCACXX"
           ["Sh.NE_Lata-007.3"]="C70LCACXX"    ["Sh.NE_Lata-033.1"]="C78YTACXX"
           ["Sh.NE_Lata-078.1"]="C78YTACXX"    ["Sh.NE_Lata-253.1"]="C78YTACXX"
           ["Sh.NE_Lata-275.2"]="C78YTACXX"    ["Sh.NE_Lata-293.1"]="C78YTACXX"
           ["Sh.NE_Lata-294.1"]="C78YTACXX"    ["Sh.NE_LibTB-009.2"]="C78YTACXX"
           ["Sh.NE_LibTB-010.1"]="C70LCACXX"   ["Sh.NE_LibTB-022.1"]="C78YTACXX"
           ["Sh.NE_LibTB-028.1"]="C78YTACXX"   ["Sh.NE_LibTB-031.1"]="C70LCACXX"
           ["Sh.NE_NG-011.1"]="C78YTACXX"      ["Sh.NE_NG-06.2"]="C78YTACXX"
           ["Sh.NE_NG-089.1"]="C70LCACXX"      ["Sh.NE_NG-236.1"]="C70LCACXX"
           ["Sh.NE_Seb-076.1"]="C78YTACXX"     ["Sh.NE_Seb-078.2"]="C78YTACXX"
           ["Sh.NE_Seb-081.2"]="C78YTACXX"     ["Sh.NE_Tiag-272.1"]="C70LCACXX"
           ["Sh.NE_YK-029.2"]="C70LCACXX"      ["Sh.NE_YK-069.1"]="C78YTACXX"
           ["Sh.NE_YK-099.2"]="C70LCACXX"      ["Sh.NE_YK-248.2"]="C78YTACXX"
           ["Sh.NE_Youri-069.2"]="C70LCACXX"   ["Sh.NE_Youri-091.3"]="C70LCACXX"
           ["Sh.TZ_PEM0063.1"]="C78YTACXX"     ["Sh.TZ_PEM0075.1"]="C78YTACXX"
           ["Sh.TZ_PEM0076.1"]="C70LCACXX"     ["Sh.TZ_PEM0079.1"]="C78YTACXX"
           ["Sh.TZ_PEM0089.2"]="C78YTACXX"     ["Sh.TZ_PEM0094.2"]="C78YTACXX"
           ["Sh.TZ_PEM0099.2"]="C78YTACXX"     ["Sh.TZ_PEM0103.1"]="C70LCACXX"
           ["Sh.TZ_PEM0104.1"]="C70LCACXX"     ["Sh.TZ_PEM0106.2"]="C78YTACXX"
           ["Sh.TZ_PEM0108.1"]="C78YTACXX"     ["Sh.TZ_PEM0110.1"]="C78YTACXX"
           ["Sh.TZ_PEM0114.3"]="C70LCACXX"     ["Sh.TZ_PEM0115.4"]="C78YTACXX"
           ["Sh.TZ_PEM0120.1"]="C78YTACXX"     ["Sh.TZ_PEM0125.1"]="C78YTACXX"
           ["Sh.TZ_PEM0126.1"]="C78YTACXX"     ["Sh.TZ_PEM0127.1"]="C78YTACXX"
           ["Sh.TZ_PEM0128.1"]="C78YTACXX"     ["Sh.TZ_PEM0130.1"]="C78YTACXX"
           ["Sh.TZ_PEM0133.1"]="C78YTACXX"     ["Sh.TZ_PEM0139.2"]="C70LCACXX"
           ["Sh.TZ_PEM0145.3"]="C70LCACXX"     ["Sh.TZ_PEM0154.1"]="C78YTACXX"
           ["Sh.TZ_PEM0157.3"]="C70LCACXX"     ["Sh.TZ_PEM0166.1"]="C78YTACXX"
           ["Sh.TZ_PEM0171.1"]="C78YTACXX"     ["Sh.TZ_UNG0006.1"]="C78YTACXX"
           ["Sh.TZ_UNG0038.1"]="C78YTACXX"     ["Sh.TZ_UNG0076.1"]="C70LCACXX"
           ["Sh.TZ_UNG0077.1"]="C70LCACXX"     ["Sh.TZ_UNG0078.1"]="C70LCACXX"
           ["Sh.TZ_UNG0087.2"]="C78YTACXX"     ["Sh.TZ_UNG0089.3"]="C78YTACXX"
           ["Sh.TZ_UNG0092.3"]="C78YTACXX"     ["Sh.TZ_UNG0099.1"]="C70LCACXX"
           ["Sh.TZ_UNG0102.1"]="C78YTACXX"     ["Sh.TZ_UNG0111.1"]="C78YTACXX"
           ["Sh.TZ_UNG0117.1"]="C78YTACXX"     ["Sh.TZ_UNG0121.1"]="C70LCACXX"
           ["Sh.TZ_UNG0125.3"]="C70LCACXX"     ["Sh.TZ_UNG0127.1"]="C78YTACXX"
           ["Sh.TZ_UNG0129.2"]="C70LCACXX"     ["Sh.TZ_UNG0134.1"]="C70LCACXX"
           ["Sh.TZ_UNG0137.3"]="C78YTACXX"     ["Sh.TZ_UNG0139.1"]="C70LCACXX"
           ["Sh.TZ_UNG0142.2"]="C70LCACXX"     ["Sh.TZ_UNG0146.1"]="C70LCACXX"
           ["ERR103051"]="1"                   ["ERR119612"]="1" 
           ["ERR119613"]="1"                   ["ERR310940"]="1"
           ["ERR539850"]="1"                   ["ERR539851"]="1" 
           ["ERR539852"]="1"                   ["ERR539853"]="1" 
           ["ERR539854"]="1"                   ["ERR539855"]="1" 
           ["ERR539856"]="1"                   ["ERR539857"]="1"        

)

RG_LANE=(   ["ERR119622"]="3"               ["ERR103048"]="4"   
            ["ERR119623"]="2"               ["Sm.BR_1278.1"]="2"
            ["Sm.BR_0447.1"]="1"            ["Sm.BR_2039.1"]="1"
            ["ERR310937"]="6"               ["ERR037800"]="6"
            ["ERR084970"]="5"               ["SRR433865"]="1"   
            ["Sh.NE_Dai-002.1"]="4"         ["Sh.NE_Dai-010.1"]="1"
            ["Sh.NE_Dai-013.3"]="2"         ["Sh.NE_Dai-031.1"]="3"
            ["Sh.NE_Dai-033.1"]="3"         ["Sh.NE_Dai-044.1"]="2"
            ["Sh.NE_Dai-045.1"]="3"         ["Sh.NE_Dai-051.1"]="4"
            ["Sh.NE_Dai-074.1"]="3"         ["Sh.NE_Dai-146.1"]="4"
            ["Sh.NE_DaiCP-233.1"]="2"       ["Sh.NE_DaiCP-276.1"]="2"
            ["Sh.NE_DaiCP-277.2"]="1"       ["Sh.NE_Doki-029.1"]="3"
            ["Sh.NE_Kar-001.1"]="4"         ["Sh.NE_Kar-002.1"]="1"
            ["Sh.NE_Kar-076.1"]="3"         ["Sh.NE_Kar-096.2"]="2"
            ["Sh.NE_Kar-241.1"]="1"         ["Sh.NE_Kar-241.2"]="3"
            ["Sh.NE_Kar-281.1"]="2"         ["Sh.NE_Kar-37.2"]="1"
            ["Sh.NE_Lata-007.3"]="2"        ["Sh.NE_Lata-033.1"]="2"
            ["Sh.NE_Lata-078.1"]="4"        ["Sh.NE_Lata-253.1"]="1"
            ["Sh.NE_Lata-275.2"]="1"        ["Sh.NE_Lata-293.1"]="3"
            ["Sh.NE_Lata-294.1"]="4"        ["Sh.NE_LibTB-009.2"]="4"
            ["Sh.NE_LibTB-010.1"]="1"       ["Sh.NE_LibTB-022.1"]="2"
            ["Sh.NE_LibTB-028.1"]="4"       ["Sh.NE_LibTB-031.1"]="1"
            ["Sh.NE_NG-011.1"]="4"          ["Sh.NE_NG-06.2"]="4"
            ["Sh.NE_NG-089.1"]="2"          ["Sh.NE_NG-236.1"]="2"
            ["Sh.NE_Seb-076.1"]="4"         ["Sh.NE_Seb-078.2"]="2"
            ["Sh.NE_Seb-081.2"]="2"         ["Sh.NE_Tiag-272.1"]="1"
            ["Sh.NE_YK-029.2"]="1"          ["Sh.NE_YK-069.1"]="2"
            ["Sh.NE_YK-099.2"]="1"          ["Sh.NE_YK-248.2"]="2"
            ["Sh.NE_Youri-069.2"]="1"       ["Sh.NE_Youri-091.3"]="2"
            ["Sh.TZ_PEM0063.1"]="1"         ["Sh.TZ_PEM0075.1"]="1"
            ["Sh.TZ_PEM0076.1"]="1"         ["Sh.TZ_PEM0079.1"]="1"
            ["Sh.TZ_PEM0089.2"]="2"         ["Sh.TZ_PEM0094.2"]="3"
            ["Sh.TZ_PEM0099.2"]="4"         ["Sh.TZ_PEM0103.1"]="2"
            ["Sh.TZ_PEM0104.1"]="1"         ["Sh.TZ_PEM0106.2"]="3"
            ["Sh.TZ_PEM0108.1"]="3"         ["Sh.TZ_PEM0110.1"]="4"
            ["Sh.TZ_PEM0114.3"]="1"         ["Sh.TZ_PEM0115.4"]="1"
            ["Sh.TZ_PEM0120.1"]="1"         ["Sh.TZ_PEM0125.1"]="2"
            ["Sh.TZ_PEM0126.1"]="2"         ["Sh.TZ_PEM0127.1"]="1"
            ["Sh.TZ_PEM0128.1"]="2"         ["Sh.TZ_PEM0130.1"]="3"
            ["Sh.TZ_PEM0133.1"]="4"         ["Sh.TZ_PEM0139.2"]="2"
            ["Sh.TZ_PEM0145.3"]="1"         ["Sh.TZ_PEM0154.1"]="3"
            ["Sh.TZ_PEM0157.3"]="2"         ["Sh.TZ_PEM0166.1"]="3"
            ["Sh.TZ_PEM0171.1"]="3"         ["Sh.TZ_UNG0006.1"]="3"
            ["Sh.TZ_UNG0038.1"]="2"         ["Sh.TZ_UNG0076.1"]="2"
            ["Sh.TZ_UNG0077.1"]="2"         ["Sh.TZ_UNG0078.1"]="1"
            ["Sh.TZ_UNG0087.2"]="1"         ["Sh.TZ_UNG0089.3"]="4"
            ["Sh.TZ_UNG0092.3"]="4"         ["Sh.TZ_UNG0099.1"]="1"
            ["Sh.TZ_UNG0102.1"]="2"         ["Sh.TZ_UNG0111.1"]="1"
            ["Sh.TZ_UNG0117.1"]="1"         ["Sh.TZ_UNG0121.1"]="1"
            ["Sh.TZ_UNG0125.3"]="2"         ["Sh.TZ_UNG0127.1"]="2"
            ["Sh.TZ_UNG0129.2"]="2"         ["Sh.TZ_UNG0134.1"]="1"
            ["Sh.TZ_UNG0137.3"]="1"         ["Sh.TZ_UNG0139.1"]="2"
            ["Sh.TZ_UNG0142.2"]="2"         ["Sh.TZ_UNG0146.1"]="1"
            ["ERR103051"]="1"                   ["ERR119612"]="1" 
            ["ERR119613"]="1"                   ["ERR310940"]="1"
            ["ERR539850"]="1"                   ["ERR539851"]="1" 
            ["ERR539852"]="1"                   ["ERR539853"]="1" 
            ["ERR539854"]="1"                   ["ERR539855"]="1" 
            ["ERR539856"]="1"                   ["ERR539857"]="1"        
        )

        

RG_INDEX=(  ["ERR119622"]="ERR119622"           ["ERR103048"]="ERR103048" 
            ["ERR119623"]="ERR119623"           ["Sm.BR_1278.1"]="CAATGGAA" 
            ["Sm.BR_0447.1"]="CTGTAGCC"         ["Sm.BR_2039.1"]="AACGCTTA"
            ["ERR310937"]="ERR310937"           ["ERR037800"]="ERR037800"
            ["ERR084970"]="ERR084970"           ["SRR433865"]="SRR433865"
            ["Sh.NE_Dai-002.1"]="GATAGACA"      ["Sh.NE_Dai-010.1"]="AACCGAGA"
            ["Sh.NE_Dai-013.3"]="AGGCTAAC"      ["Sh.NE_Dai-031.1"]="CCGACAAC"
            ["Sh.NE_Dai-033.1"]="GAGTTAGC"      ["Sh.NE_Dai-044.1"]="AGCCATGC"
            ["Sh.NE_Dai-045.1"]="GATGAATC"      ["Sh.NE_Dai-051.1"]="TAGGATGA"
            ["Sh.NE_Dai-074.1"]="ATTGGCTC"      ["Sh.NE_Dai-146.1"]="GCGAGTAA"
            ["Sh.NE_DaiCP-233.1"]="ATCGCGAC"    ["Sh.NE_DaiCP-276.1"]="AGCAGGAA"
            ["Sh.NE_DaiCP-277.2"]="AGTACAAG"    ["Sh.NE_Doki-029.1"]="CAAGGAGC"
            ["Sh.NE_Kar-001.1"]="GCTAACGA"      ["Sh.NE_Kar-002.1"]="ATGCCTAA"
            ["Sh.NE_Kar-076.1"]="CGGATTGC"      ["Sh.NE_Kar-096.2"]="AGTCACTA"
            ["Sh.NE_Kar-241.1"]="AACGCTTA"      ["Sh.NE_Kar-241.2"]="CTAAGGTC"
            ["Sh.NE_Kar-281.1"]="ACACAGAA"      ["Sh.NE_Kar-37.2"]="CCTCCTGA"
            ["Sh.NE_Lata-007.3"]="AAATCACC"     ["Sh.NE_Lata-033.1"]="ACGTATCA"
            ["Sh.NE_Lata-078.1"]="GTCGTAGA"     ["Sh.NE_Lata-253.1"]="AGTGGTCA"
            ["Sh.NE_Lata-275.2"]="ACCACTGT"     ["Sh.NE_Lata-293.1"]="CCATCCTC"
            ["Sh.NE_Lata-294.1"]="GCCACATA"     ["Sh.NE_LibTB-009.2"]="TCCGTCTA"
            ["Sh.NE_LibTB-010.1"]="GAGNTGAA"    ["Sh.NE_LibTB-022.1"]="CAATGGAA"
            ["Sh.NE_LibTB-028.1"]="GGTGCGAA"    ["Sh.NE_LibTB-031.1"]="CTCAATGA"
            ["Sh.NE_NG-011.1"]="GTACGCAA"       ["Sh.NE_NG-06.2"]="GTGTTCTA"
            ["Sh.NE_NG-089.1"]="AGCACCTC"       ["Sh.NE_NG-236.1"]="AATGTTGC"
            ["Sh.NE_Seb-076.1"]="GTCTGTCA"      ["Sh.NE_Seb-078.2"]="ATTGAGGA"
            ["Sh.NE_Seb-081.2"]="ATCCTGTA"      ["Sh.NE_Tiag-272.1"]="CGANCTTA"
            ["Sh.NE_YK-029.2"]="CAGNGTTA"       ["Sh.NE_YK-069.1"]="ACGCTCGA"
            ["Sh.NE_YK-099.2"]="CCAGTTCA"       ["Sh.NE_YK-248.2"]="ACTATGCA"
            ["Sh.NE_Youri-069.2"]="CTGGCATA"    ["Sh.NE_Youri-091.3"]="ACCCGACC"
            ["Sh.TZ_PEM0063.1"]="CAGATCTG"      ["Sh.TZ_PEM0075.1"]="AAACATCG"
            ["Sh.TZ_PEM0076.1"]="CGCNTACA"      ["Sh.TZ_PEM0079.1"]="ACAAGCTA"
            ["Sh.TZ_PEM0089.2"]="CAACCACA"      ["Sh.TZ_PEM0094.2"]="CACCTTAC"
            ["Sh.TZ_PEM0099.2"]="TCTTCACA"      ["Sh.TZ_PEM0103.1"]="TTNACGCA"
            ["Sh.TZ_PEM0104.1"]="CTGAGCCA"      ["Sh.TZ_PEM0106.2"]="GAACAGGC"
            ["Sh.TZ_PEM0108.1"]="CCTAATCC"      ["Sh.TZ_PEM0110.1"]="GGAGAACA"
            ["Sh.TZ_PEM0114.3"]="GACTAGTA"      ["Sh.TZ_PEM0115.4"]="AACGTGAT"
            ["Sh.TZ_PEM0120.1"]="CGCTGATC"      ["Sh.TZ_PEM0125.1"]="ACAGCAGA"
            ["Sh.TZ_PEM0126.1"]="AGATCGCA"      ["Sh.TZ_PEM0127.1"]="ACATTGGC"
            ["Sh.TZ_PEM0128.1"]="AAGGTACA"      ["Sh.TZ_PEM0130.1"]="CCTCTATC"
            ["Sh.TZ_PEM0133.1"]="TATCAGCA"      ["Sh.TZ_PEM0139.2"]="TGNAACAA"
            ["Sh.TZ_PEM0145.3"]="CACTTCGA"      ["Sh.TZ_PEM0154.1"]="GACAGTGC"
            ["Sh.TZ_PEM0157.3"]="AANAGATC"      ["Sh.TZ_PEM0166.1"]="CGACACAC"
            ["Sh.TZ_PEM0171.1"]="GCCAAGAC"      ["Sh.TZ_UNG0006.1"]="ATCATTCC"
            ["Sh.TZ_UNG0038.1"]="AGAGTCAA"      ["Sh.TZ_UNG0076.1"]="TGNTGGTA"
            ["Sh.TZ_UNG0077.1"]="AGNTGTAC"      ["Sh.TZ_UNG0078.1"]="GAATCTGA"
            ["Sh.TZ_UNG0087.2"]="CTGTAGCC"      ["Sh.TZ_UNG0089.3"]="TGAAGAGA"
            ["Sh.TZ_UNG0092.3"]="GCTCGGTA"      ["Sh.TZ_UNG0099.1"]="CGACTGGA"
            ["Sh.TZ_UNG0102.1"]="CAAGACTA"      ["Sh.TZ_UNG0111.1"]="AACAACCA"
            ["Sh.TZ_UNG0117.1"]="AAGACGGA"      ["Sh.TZ_UNG0121.1"]="CCGTGAGA"
            ["Sh.TZ_UNG0125.3"]="AANCCGTC"      ["Sh.TZ_UNG0127.1"]="ACCTCCAA"
            ["Sh.TZ_UNG0129.2"]="TGGCTTCA"      ["Sh.TZ_UNG0134.1"]="CCGAGGTA"
            ["Sh.TZ_UNG0137.3"]="CATCAAGT"      ["Sh.TZ_UNG0139.1"]="AAAGACAC"
            ["Sh.TZ_UNG0142.2"]="ACAGATTC"      ["Sh.TZ_UNG0146.1"]="CATACCAA"
            ["ERR103051"]="1"                   ["ERR119612"]="1" 
            ["ERR119613"]="1"                   ["ERR310940"]="1"
            ["ERR539850"]="1"                   ["ERR539851"]="1" 
            ["ERR539852"]="1"                   ["ERR539853"]="1" 
            ["ERR539854"]="1"                   ["ERR539855"]="1" 
            ["ERR539856"]="1"                   ["ERR539857"]="1"        
        )


GENOMES=(   
    ["Sm.BR_0447.1"]="$MAN_GENOME"          ["Sm.BR_1278.1"]="$MAN_GENOME"
    ["Sm.BR_2039.1"]="$MAN_GENOME"          ["ERR103048"]="$BOV_GENOME"
    ["ERR119622"]="$BOV_GENOME"             ["ERR119623"]="$CUR_GENOME"
    ["ERR310937"]="$CUR_GENOME"             ["SRR433865"]="$HAE_GENOME"
    ["ERR037800"]="$HAE_GENOME"             ["ERR084970"]="$HAE_GENOME"
    ["Sh.NE_Dai-002.1"]="$HAE_GENOME"       ["Sh.NE_Dai-010.1"]="$HAE_GENOME"
    ["Sh.NE_Dai-013.3"]="$HAE_GENOME"       ["Sh.NE_Dai-031.1"]="$HAE_GENOME"
    ["Sh.NE_Dai-033.1"]="$HAE_GENOME"       ["Sh.NE_Dai-044.1"]="$HAE_GENOME"
    ["Sh.NE_Dai-045.1"]="$HAE_GENOME"       ["Sh.NE_Dai-051.1"]="$HAE_GENOME"
    ["Sh.NE_Dai-074.1"]="$HAE_GENOME"       ["Sh.NE_Dai-146.1"]="$HAE_GENOME"
    ["Sh.NE_DaiCP-233.1"]="$HAE_GENOME"     ["Sh.NE_DaiCP-276.1"]="$HAE_GENOME"
    ["Sh.NE_DaiCP-277.2"]="$HAE_GENOME"     ["Sh.NE_Doki-029.1"]="$HAE_GENOME"
    ["Sh.NE_Kar-001.1"]="$HAE_GENOME"       ["Sh.NE_Kar-002.1"]="$HAE_GENOME"
    ["Sh.NE_Kar-076.1"]="$HAE_GENOME"       ["Sh.NE_Kar-096.2"]="$HAE_GENOME"
    ["Sh.NE_Kar-241.1"]="$HAE_GENOME"       ["Sh.NE_Kar-241.2"]="$HAE_GENOME"
    ["Sh.NE_Kar-281.1"]="$HAE_GENOME"       ["Sh.NE_Kar-37.2"]="$HAE_GENOME"
    ["Sh.NE_Lata-007.3"]="$HAE_GENOME"      ["Sh.NE_Lata-033.1"]="$HAE_GENOME"
    ["Sh.NE_Lata-078.1"]="$HAE_GENOME"      ["Sh.NE_Lata-253.1"]="$HAE_GENOME"
    ["Sh.NE_Lata-275.2"]="$HAE_GENOME"      ["Sh.NE_Lata-293.1"]="$HAE_GENOME"
    ["Sh.NE_Lata-294.1"]="$HAE_GENOME"      ["Sh.NE_LibTB-009.2"]="$HAE_GENOME"
    ["Sh.NE_LibTB-010.1"]="$HAE_GENOME"     ["Sh.NE_LibTB-022.1"]="$HAE_GENOME"
    ["Sh.NE_LibTB-028.1"]="$HAE_GENOME"     ["Sh.NE_LibTB-031.1"]="$HAE_GENOME"
    ["Sh.NE_NG-011.1"]="$HAE_GENOME"        ["Sh.NE_NG-06.2"]="$HAE_GENOME"
    ["Sh.NE_NG-089.1"]="$HAE_GENOME"        ["Sh.NE_NG-236.1"]="$HAE_GENOME"
    ["Sh.NE_Seb-076.1"]="$HAE_GENOME"       ["Sh.NE_Seb-078.2"]="$HAE_GENOME"
    ["Sh.NE_Seb-081.2"]="$HAE_GENOME"       ["Sh.NE_Tiag-272.1"]="$HAE_GENOME"
    ["Sh.NE_YK-029.2"]="$HAE_GENOME"        ["Sh.NE_YK-069.1"]="$HAE_GENOME"
    ["Sh.NE_YK-099.2"]="$HAE_GENOME"        ["Sh.NE_YK-248.2"]="$HAE_GENOME"
    ["Sh.NE_Youri-069.2"]="$HAE_GENOME"     ["Sh.NE_Youri-091.3"]="$HAE_GENOME"
    ["Sh.TZ_PEM0063.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0075.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0076.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0079.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0089.2"]="$HAE_GENOME"       ["Sh.TZ_PEM0094.2"]="$HAE_GENOME"
    ["Sh.TZ_PEM0099.2"]="$HAE_GENOME"       ["Sh.TZ_PEM0103.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0104.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0106.2"]="$HAE_GENOME"
    ["Sh.TZ_PEM0108.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0110.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0114.3"]="$HAE_GENOME"       ["Sh.TZ_PEM0115.4"]="$HAE_GENOME"
    ["Sh.TZ_PEM0120.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0125.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0126.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0127.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0128.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0130.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0133.1"]="$HAE_GENOME"       ["Sh.TZ_PEM0139.2"]="$HAE_GENOME"
    ["Sh.TZ_PEM0145.3"]="$HAE_GENOME"       ["Sh.TZ_PEM0154.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0157.3"]="$HAE_GENOME"       ["Sh.TZ_PEM0166.1"]="$HAE_GENOME"
    ["Sh.TZ_PEM0171.1"]="$HAE_GENOME"       ["Sh.TZ_UNG0006.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0038.1"]="$HAE_GENOME"       ["Sh.TZ_UNG0076.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0077.1"]="$HAE_GENOME"       ["Sh.TZ_UNG0078.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0087.2"]="$HAE_GENOME"       ["Sh.TZ_UNG0089.3"]="$HAE_GENOME"
    ["Sh.TZ_UNG0092.3"]="$HAE_GENOME"       ["Sh.TZ_UNG0099.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0102.1"]="$HAE_GENOME"       ["Sh.TZ_UNG0111.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0117.1"]="$HAE_GENOME"       ["Sh.TZ_UNG0121.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0125.3"]="$HAE_GENOME"       ["Sh.TZ_UNG0127.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0129.2"]="$HAE_GENOME"       ["Sh.TZ_UNG0134.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0137.3"]="$HAE_GENOME"       ["Sh.TZ_UNG0139.1"]="$HAE_GENOME"
    ["Sh.TZ_UNG0142.2"]="$HAE_GENOME"       ["Sh.TZ_UNG0146.1"]="$HAE_GENOME"
    )




