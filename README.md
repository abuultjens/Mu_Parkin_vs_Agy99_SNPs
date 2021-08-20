# Mu_Parkin_vs_Agy99_SNPs

Want to show that the SNPs detected when using Parkin are better than AGY99

Using 478 isolates with the updated epi labels

### WD

    # Agy99_chr-p
    /home/buultjensa/2020_Mu/snippy_v4.6.0/Agy99_chr-p
    
    fa -v 478-Mu_s4.6.0_Agy99-chr-p.aln
    (stdin)                   no=479 bp=2212980 ok=2212980 Ns=0 gaps=0 min=4620 avg=4620 max=4620 N50=4620
    
    fa -v 478-Mu_s4.6.0_Agy99-chr-p.noref.aln
    (stdin)                   no=478 bp=249038 ok=249038 Ns=0 gaps=0 min=521 avg=521 max=521 N50=521

    # Mu_Parkin_2021_chr-p
    /home/buultjensa/2020_Mu/snippy_v4.6.0/Mu_Parkin_2021_chr-p

    fa -v 478-Mu_s4.6.0_Parkin-chr-p.aln
    (stdin)                   no=479 bp=150406 ok=150406 Ns=0 gaps=0 min=314 avg=314 max=314 N50=314

    fa -v 478-Mu_s4.6.0_Parkin-chr-p.noref.aln
    (stdin)                   no=478 bp=136708 ok=136708 Ns=0 gaps=0 min=286 avg=286 max=286 N50=286
    
### Extracting 100bp up and down of SNP positions

#### PARKIN

    # count number of SNPs in Parkin file
    wc -l 478_Parkin_chr_p_pos.txt 
    285
    
    ### extract parkin SNP regions to file
    # make seq file with chr on single line
    samtools faidx M_ulcerans_JKD8049.fa
    samtools faidx M_ulcerans_JKD8049.fa Mulcerans_JKD8049 | grep -v ">" | tr -d '\n' > Mulcerans_JKD8049_1LINE.seq
    echo "" >> Mulcerans_JKD8049_1LINE.seq
    for TAXA in $(cat 478_Parkin_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}         
        echo ">SNP-${TAXA}_Parkin-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS.fa    
        cut -b ${LOW}-${HIGH} Mulcerans_JKD8049_1LINE.seq >> Parkin_Agy99_SNP-REGIONS.fa
    done

#### AGY99

    # count number of SNPs in Agy99 file
    wc -l 478_Agy99_chr_p_pos.txt 
    520
    
    ### extract Agy99 SNP regions to file
    # make seq file with chr on single line
    samtools faidx Agy99-chr-p.fa
    samtools faidx Agy99-chr-p.fa CP000325.1 | grep -v ">" | tr -d '\n' > Agy99-chr-p_1LINE.seq
    echo "" >> Agy99-chr-p_1LINE.seq
    for TAXA in $(cat 478_Agy99_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}          
        echo ">SNP-${TAXA}_Agy99-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS.fa
        cut -b ${LOW}-${HIGH} Agy99-chr-p_1LINE.seq >> Parkin_Agy99_SNP-REGIONS.fa      
    done 
         
### clustering SNP regions with cd-hit-est

    cd-hit-est -i Parkin_Agy99_SNPs.fa -o Parkin_Agy99_SNP-REGIONS_cd-hit-est_c-0.8 -d 120 -c 0.8
    
518 clusters:  
230 Agy99 singletons  
11 Parkin singletons  
271 clusters with two members  
2 clusters with three members  
4 clusters with four members  

### SNP 25580 in Agy99

The 201 bp region containing SNP 25580 in Agy99 (100 bases downstream and 100 bases upstream) did not have a homolog among the Parkin SNP regions at clustering ID of 0.8 and was one of the Agy99 singletons. For this site Agy99 has a T and in the core.tab file 475 isolates have a C and 3 isolates have a T. I manually checked the bam files for the three isolates with the C and three with the T allele reported in the core.tab file.

#### 2020-12845 has a T in core.tab  
-bam file shows a C SNP when mapped to Agy99  
-40x cov for C    

    # grepping for site 25580 in snps.vcf
    grep "25580" 2020-12845/snps.vcf
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2020-12845/snps.csv
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2020-12845/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

#### 2015-104 has a T in core.tab  
-bam file shows a C SNP when mapped to Agy99  
-54x cov for C  

    # grepping for site 25580 in snps.vcf
    grep "25580" 2015-104/snps.vcf
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2015-104/snps.csv
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2015-104/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

#### 2018-11426 has a T in core.tab  
-bam file shows a C SNP when mapped to Agy99  
-85x cov for C  

    # grepping for site 25580 in snps.vcf
    grep "25580" 2018-11426/snps.vcf
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2018-11426/snps.csv
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2018-11426/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

#### 2020-12844 has a C in core.tab  
-bam file shows a C SNP when mapped to Agy99  
-34x cov for C

    # grepping for site 25580 in snps.vcf
    grep "25580" 2020-12844/snps.vcf
    CP000325.1      25580   .       T       C       754.801 .       AB=0;AO=29;DP=31;QA=911;QR=0;RO=0;TYPE=snp      GT:DP:RO:QR:AO:QA:GL        1/1:31:0:0:29:911:-81.968,-8.72987,0
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2020-12844/snps.csv
    CP000325.1,25580,snp,T,C,C:29 T:0
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2020-12844/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C
    
#### 16156690
-bam file shows a C SNP when mapped to Agy99  
-101x cov for C  

    # grepping for site 25580 in snps.vcf
    grep "25580" 16156690/snps.vcf
    CP000325.1      25580   .       T       C       2542.76 .       AB=0;AO=97;DP=101;QA=2869;QR=0;RO=0;TYPE=snp    GT:DP:RO:QR:AO:QA:GL       1/1:101:0:0:97:2869:-258.236,-29.1999,0
    
    # grepping for site 25580 in snps.csv
    grep "25580" 16156690/snps.csv
    CP000325.1,25580,snp,T,C,C:97 T:0
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 16156690/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

#### 17104468
-bam file shows a C SNP when mapped to Agy99  
-85x cov for C  

    # grepping for site 25580 in snps.vcf
    grep "25580" 17104468/snps.vcf
    CP000325.1      25580   .       T       C       2232.85 .       AB=0;AO=85;DP=85;QA=2517;QR=0;RO=0;TYPE=snp     GT:DP:RO:QR:AO:QA:GL       1/1:85:0:0:85:2517:-226.748,-25.5876,0

    # grepping for site 25580 in snps.csv
    grep "25580" 17104468/snps.csv
    CP000325.1,25580,snp,T,C,C:85 T:0
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 17104468/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

### using samtools tview to check the consensus base directly from the bam file

I put the above samtools tview command in a loop and checked all SNP positions for all isolates and used Excel to determine that 237 sites became invariant, leaving 283 actual SNPs.

### Extracting 100bp up and down of SNP positions from the samtools tview verified Agy99 alignment

#### PARKIN

    # count number of SNPs in Parkin file
    wc -l 478_Parkin_chr_p_pos.txt 
    285
    
    # extract parkin SNP regions to file
    for TAXA in $(cat 478_Parkin_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}         
        echo ">SNP-${TAXA}_Parkin-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa    
        cut -b ${LOW}-${HIGH} ../../Mu_Parkin_2021_chr-p/Mulcerans_JKD8049_1LINE.seq >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa
    done

#### AGY99

    # count number of SNPs in Agy99 file
    wc -l 478_Agy99_chr_p_pos_AFTER-TVIEW.txt 
    283
    
    # extract Agy99 SNP regions to file
    for TAXA in $(cat 478_Agy99_chr_p_pos_AFTER-TVIEW.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}          
        echo ">SNP-${TAXA}_Agy99-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa
        cut -b ${LOW}-${HIGH} ../Agy99-chr-p_1LINE.seq >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa      
    done 
         
### clustering SNP regions with cd-hit-est

    cd-hit-est -i Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa   -o Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW_cd-hit-est_c-0.8 -d 120 -c 0.8
    
292 clusters:  
10 Agy99 singletons  
12 Parkin singletons  
267 clusters with two members  
3 clusters with four members  

The SNP region clustering provides a way to align SNPs from differnet references. Depending on how distant the references are there will be some SNPs that are specific to each ref. It was the high number of SNPs that were specific to Agy99 (the Agy99 singleton clusters) that alerted me to investigate them and then find that they are false SNPs.


### stats

Agy99 verified SNP alignment (283 SNPs):
197 were singletons SNPs, the minor allel only occured in a single isolate in the set
86 had a minor allele that occured in two or more isolates

Parkin SNP alignment (285 SNPs):
200 were singletons, the minor allel only occured in a single isolate in the set
85 had a minor allele that occured in two or more isolates



