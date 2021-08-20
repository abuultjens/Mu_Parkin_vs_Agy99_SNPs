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
    
    # extract parkin SNP regions to file
    for TAXA in $(cat 478_Parkin_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}         
        echo ">SNP-${TAXA}_Parkin-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS.fa    
        cut -b ${LOW}-${HIGH} ../../Mu_Parkin_2021_chr-p/Mulcerans_JKD8049_1LINE.seq >> Parkin_Agy99_SNP-REGIONS.fa
    done

#### AGY99

    # count number of SNPs in Agy99 file
    wc -l 478_Agy99_chr_p_pos.txt 
    520
    
    # extract Agy99 SNP regions to file
    for TAXA in $(cat 478_Agy99_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}          
        echo ">SNP-${TAXA}_Agy99-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS.fa
        cut -b ${LOW}-${HIGH} ../Agy99-chr-p_1LINE.seq >> Parkin_Agy99_SNP-REGIONS.fa      
    done 
         
### clustering SNP regions with cd-hit-est

    cd-hit-est -i Parkin_Agy99_SNPs.fa -o Parkin_Agy99_SNP-REGIONS_cd-hit-est_c-0.8 -d 120 -c 0.8
    
518 clusters:  
230 Agy99 singletons  
11 Parkin singletons  
271 clusters with two members  
2 clusters with two members  
4 clusters with four members  

    

    # here is a SNP region in Agy99 that had no homology in Parkin
    SNP-25580_Agy99-chr-p_201.fa 
    GGTAGCCCTGGTCGTAGCCGCCACCCTGGTCGGGGTAGCCCGGGCGCTGCGCCTCAGGGGGTGCCGACGGAGCCGCGGGCGCTGAATAGGCGCCCTCGTCTTGGCGTGCGGGCTCGCGGCCGTATTCGCCGTAGCCCCCGTATCCGGGCTGACCGCCCGGCGGCTGGCCGTAGCCGCCACCCTGGCGGTAACCCTGGTCGT

I blasted this region and found a match in Parkin chr:25479-25679  

2020-12845 has a T in core.tab  
-looks like a C when mapped to Parkin  
-has a C SNP when mapped to Agy99  

2020-12844 has a C in core.tab  
-looks like a C when mapped to Parkin  
-Has a C SNP when mapped to Agy99  

2015-104 has a T in core.tab  
-Has a C SNP when mapped to Agy99  

2018-11426 has a T in core.tab  
-Has a C SNP when mapped to Agy99  

### using samtools tview to check SNPs

    for TAXA in $(cat $1); do
        # write header line
        echo "INDEX,POSITION" > cov-checker_Agy99-chr-p_4564243.csv
        # get base
        BASE=`sh ~/shell_scripts/Snippy/cov-checker.sh ../Agy99-chr-p.fa CP000325.1 4564243 ../${TAXA}/snps.bam | tail -2 | head -1 | cut -f 2 -d ' ' | cut -f 2 -d '[' | cut -f 1 -d ']' | cut -f 1 -d ','`
        # write body line
        echo "${TAXA},${BASE}" >> cov-checker_Agy99-chr-p_4564243.csv
    done

### Used Excel to determine invariant positions in alignment

Agy99_chr_p 478 had 283 SNPs

### Extracting 100bp up and down of SNP positions

#### PARKIN

    # count number of SNPs in Parkin file
    wc -l 478_Parkin_chr_p_pos.txt 
    285
    
    # extract parkin SNP regions to file
    for TAXA in $(cat 478_Parkin_chr_p_pos.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}         
        echo ">SNP-${TAXA}_Parkin-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS_AFTER-COV-CHECKER.fa    
        cut -b ${LOW}-${HIGH} ../../Mu_Parkin_2021_chr-p/Mulcerans_JKD8049_1LINE.seq >> Parkin_Agy99_SNP-REGIONS_AFTER-COV-CHECKER.fa
    done

#### AGY99

    # count number of SNPs in Agy99 file
    wc -l 478_Agy99_chr_p_pos_AFTER-COV-CHECKER.txt 
    283
    
    # extract Agy99 SNP regions to file
    for TAXA in $(cat 478_Agy99_chr_p_pos_AFTER-COV-CHECKER.txt); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}          
        echo ">SNP-${TAXA}_Agy99-chr-p_201_${LOW}-${HIGH}" >> Parkin_Agy99_SNP-REGIONS_AFTER-COV-CHECKER.fa
        cut -b ${LOW}-${HIGH} ../Agy99-chr-p_1LINE.seq >> Parkin_Agy99_SNP-REGIONS_AFTER-COV-CHECKER.fa      
    done 
         
### clustering SNP regions with cd-hit-est

    cd-hit-est -i Parkin_Agy99_SNP-REGIONS_AFTER-COV-CHECKER.fa   -o Parkin_Agy99_SNP-REGIONS_AFTER-COV-CHECKER_cd-hit-est_c-0.8 -d 120 -c 0.8
    
292 clusters:  
10 Agy99 singletons  
12 Parkin singletons  
267 clusters with two members  
3 clusters with four members  

There about five Parkin and Agy99 singletons

The SNP region clustering provides a way to align SNPs from differnet references. Depending on how distant the references are there will be some SNPs that are specific to each ref. It was the high number of SNPs that were specific to Agy99 (the Agy99 singleton clusters) that alerted me to investigate them and then find that they are false SNPs.


