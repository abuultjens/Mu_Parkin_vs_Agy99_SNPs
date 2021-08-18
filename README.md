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

    for TAXA in $(cat $1); do
        NUMB=100
        let LOW=${TAXA}-${NUMB}
        let HIGH=${TAXA}+${NUMB}        
        echo ">SNP-${TAXA}_Agy99-chr-p_201_${LOW}-${HIGH}" > Parkin_Agy99_SNPs.fa     
        cut -b ${LOW}-${HIGH} Agy99-chr-p_1LINE.seq >> Parkin_Agy99_SNPs.fa
    done 
    
### clustering of regions with cd-hit-est
    cd-hit-est -i Parkin_Agy99_SNPs.fa -o out_cd-hit-est_c-0.8 -d 120 -c 0.8
    
    
    # here is a SNP region in Agy99 that had no homology in Parkin
    SNP-25580_Agy99-chr-p_201.fa 
    GGTAGCCCTGGTCGTAGCCGCCACCCTGGTCGGGGTAGCCCGGGCGCTGCGCCTCAGGGGGTGCCGACGGAGCCGCGGGCGCTGAATAGGCGCCCTCGTCTTGGCGTGCGGGCTCGCGGCCGTATTCGCCGTAGCCCCCGTATCCGGGCTGACCGCCCGGCGGCTGGCCGTAGCCGCCACCCTGGCGGTAACCCTGGTCGT

I blasted this region and found a match in Parkin chr:25479-25679
2020-12845 has a SNP (T)





