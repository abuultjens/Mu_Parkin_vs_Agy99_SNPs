# Mu_Parkin_vs_Agy99_SNPs

I have found 237 false SNPs in snippy and snippy-core files when mapping 478 Vic Mu isolates to the Agy99_chr_p ref. I did not find any false SNPs when mapping the same isolate set against the Parkin_chr_p ref. Snippy v4.6.0 was used. I used samtools tview to verify each SNP for each isolate directly from the bam files. The tview verified SNP alignments have approximatly the same number of sites, meeting expectations for this pathogen.  

### WD  

    # Agy99_chr_p
    /home/buultjensa/2020_Mu/snippy_v4.6.0/Agy99_chr-p
    
    # Parkin_chr_p
    /home/buultjensa/2020_Mu/snippy_v4.6.0/Mu_Parkin_2021_chr-p

### SNP counts
    
Mapping a set of 478 Mu Victorian isolates with snippv v4.6.0 against Agy99_chr_p gave 520 SNPs while mapping the same set to Parkin_chr_p gave 285 SNPs. My expectation was that I should get more SNPs with Parkin as the ref as it is a Vic isolate and Agy99 is around 5,000 SNPs divergent to the Vic Mu population. The closer related Vic ref should have more potential for SNP sites, as more reads from the Vic isolate set should align due to it having regions of difference in common compared to the distant Agy99.
    
    # 478 Mu mapped to Agy99_chr_p 
    520 SNPs (excluding a single SNP on the plasmid)
    
    # 478 Mu mapped to Parkin_chr_p
    285 SNPs (excluding a single SNP on the plasmid)

### Extracting and clustering regions spanning 100bp up and down of SNP positions

I wanted to inspect what SNPs are the same between the references and determine those that are specific to the Agy99 alignment. This was difficult to do as they are from different references. In order to compare the SNPs between these references I extracted regions containing the SNP sites to provide flanking context for a clustering comparison. Here I took 100 bp up and downstream of the SNP site, totalling 201 bp per SNP region. I then used cd-hit-est to cluster the SNP regions at ID cutoff of 0.8.

    # WD
    /home/buultjensa/2020_Mu/snippy_v4.6.0/Agy99_chr-p/Parkin_vs_Agy99

#### Parkin

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

#### Agy99

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
         
#### Clustering SNP regions with cd-hit-est

    cd-hit-est -i Parkin_Agy99_SNPs.fa -o Parkin_Agy99_SNP-REGIONS_cd-hit-est_c-0.8 -d 120 -c 0.8
    
**518 clusters:**  
230 Agy99 singletons  
11 Parkin singletons  
271 clusters with two members  
2 clusters with three members  
4 clusters with four members  

My expectation was that majority of clusters should contain one SNP region from both Agy99 and parkin (close to a 1:1 merge of SNPs). There was a surprisingly high number of singleton clusters that had just a single Agy99 SNP region.

### Investigating the Agy99 singleton SNP region clusters: SNP 25580

The 201 bp region containing SNP 25580 in Agy99 (100 bases downstream and 100 bases upstream) did not have a homolog among the Parkin SNP regions at clustering ID of 0.8 and was one of the 230 Agy99 singletons. For this site the Agy99 ref had a T and in the core.tab file 475 isolates had a C and 3 isolates had a T. I manually checked the bam files for the three isolates with the C and three with the T allele, as reported in the core.tab file.

    # WD
    /home/buultjensa/2020_Mu/snippy_v4.6.0/Agy99_chr-p

**2020-12845 has a T in core.tab for site 25580**  
-bam file shows a C SNP when mapped to Agy99  
-40x cov for C (samtools tview)    

    # grepping for site 25580 in snps.vcf
    grep "25580" 2020-12845/snps.vcf
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2020-12845/snps.csv
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2020-12845/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

**2015-104 has a T in core.tab for site 25580**  
-bam file shows a C SNP when mapped to Agy99  
-54x cov for C (samtools tview)  

    # grepping for site 25580 in snps.vcf
    grep "25580" 2015-104/snps.vcf
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2015-104/snps.csv
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2015-104/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

**2018-11426 has a T in core.tab for site 25580**  
-bam file shows a C SNP when mapped to Agy99  
-85x cov for C (samtools tview)  

    # grepping for site 25580 in snps.vcf
    grep "25580" 2018-11426/snps.vcf
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2018-11426/snps.csv
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2018-11426/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

**2020-12844 has a C in core.tab for site 25580**  
-bam file shows a C SNP when mapped to Agy99  
-34x cov for C (samtools tview)

    # grepping for site 25580 in snps.vcf
    grep "25580" 2020-12844/snps.vcf
    CP000325.1      25580   .       T       C       754.801 .       AB=0;AO=29;DP=31;QA=911;QR=0;RO=0;TYPE=snp      GT:DP:RO:QR:AO:QA:GL        1/1:31:0:0:29:911:-81.968,-8.72987,0
    
    # grepping for site 25580 in snps.csv
    grep "25580" 2020-12844/snps.csv
    CP000325.1,25580,snp,T,C,C:29 T:0
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 2020-12844/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C
    
**16156690 has a C in core.tab for site 25580**  
-bam file shows a C SNP when mapped to Agy99  
-101x cov for C (samtools tview)  

    # grepping for site 25580 in snps.vcf
    grep "25580" 16156690/snps.vcf
    CP000325.1      25580   .       T       C       2542.76 .       AB=0;AO=97;DP=101;QA=2869;QR=0;RO=0;TYPE=snp    GT:DP:RO:QR:AO:QA:GL       1/1:101:0:0:97:2869:-258.236,-29.1999,0
    
    # grepping for site 25580 in snps.csv
    grep "25580" 16156690/snps.csv
    CP000325.1,25580,snp,T,C,C:97 T:0
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 16156690/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

**17104468 has a C in core.tab for site 25580**  
-bam file shows a C SNP when mapped to Agy99  
-85x cov for C (samtools tview)  

    # grepping for site 25580 in snps.vcf
    grep "25580" 17104468/snps.vcf
    CP000325.1      25580   .       T       C       2232.85 .       AB=0;AO=85;DP=85;QA=2517;QR=0;RO=0;TYPE=snp     GT:DP:RO:QR:AO:QA:GL       1/1:85:0:0:85:2517:-226.748,-25.5876,0

    # grepping for site 25580 in snps.csv
    grep "25580" 17104468/snps.csv
    CP000325.1,25580,snp,T,C,C:85 T:0
    
    # using samtools tview to check the consensus base directly from the bam file
    samtools tview -d T 17104468/snps.bam Agy99-chr-p.fa -p CP000325.1:25580 | cut -b 1 | head -3 | tail -1
    C

### Using samtools tview to check the consensus base directly from the bam file

I ran the same samtools tview command in a loop and checked all SNP positions for all isolates and used Excel to determine that 237 sites became invariant (were false SNPs), leaving 283 actual SNPs. I ran the same tview script on the parkin SNP alignment and no sites were found to become invariant.

### Extracting 100bp up and down of SNP positions from the samtools tview verified Agy99 alignment

#### Parkin

    # count number of SNPs in Parkin file
    wc -l 478_Parkin_chr_p_pos.txt 
    285
    
    # extract parkin SNP regions to file
    for TAXA in $(cat 478_Parkin_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}         
        echo ">SNP-${TAXA}_Parkin-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa    
        cut -b ${LOW}-${HIGH} Mulcerans_JKD8049_1LINE.seq >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa
    done

#### Agy99

    # count number of SNPs in Agy99 file
    wc -l 478_Agy99_chr_p_pos_AFTER-TVIEW.txt 
    283
    
    # extract Agy99 SNP regions to file
    for TAXA in $(cat 478_Agy99_chr_p_pos_AFTER-TVIEW.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}          
        echo ">SNP-${TAXA}_Agy99-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa
        cut -b ${LOW}-${HIGH} Agy99-chr-p_1LINE.seq >> Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa      
    done 
         
### Clustering tview verified SNP regions with cd-hit-est

    cd-hit-est -i Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW.fa -o Parkin_Agy99_SNP-REGIONS_AFTER-TVIEW_cd-hit-est_c-0.8 -d 120 -c 0.8
    
**292 clusters:**  
10 Agy99 singletons  
12 Parkin singletons  
267 clusters with two members  
3 clusters with four members  

This SNP region clustering approach provides a way to align SNPs from different references. Depending on how distant the references are there will be some SNPs that are specific to each ref (singleton clusters). It was the high number of SNPs that were specific to Agy99 (the Agy99 singleton clusters) that alerted me to investigate them and find out that they were false. After removing the false SNPs specific to Agy99, the majority of SNP regions grouped into clusters with two members, one from each ref (1:1 match between SNPs from each ref), which is what is to be expected for most of the SNPs. 

### SNP stats

**Agy99 verified SNP alignment (283 SNPs):**  
197 were singletons SNPs, the minor allele only occurred in a single isolate in the set  
86 had a minor allele that occurred in two or more isolates  

**Parkin SNP alignment (285 SNPs):**  
200 were singletons, the minor allele only occurred in a single isolate in the set  
85 had a minor allele that occurred in two or more isolates  

### Agy99 singletons after tview

Agy99 SNP 398430

homologus region in Parkin is: 467346-467146

***2020-12842 has a G in tview-core.tab***  



***2020-12844 has a A in tview-core.tab***  




