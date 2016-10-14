# SNP_INDEL_Annotations

Annotates coding SNPs and INDELs from a genome sequence (fasta) and gene annotation file (gtf).
Currently only works for Ensemble annotation files. However, it would not be hard to adjust for something like refSeq.
Pretty slow unfortunately, but it's accurate at least and takes into account snps within the same codon.
This was written a couple of years ago, but I did end up comparing the results to snpEff eventually and they were the same,
except snpEff did not properly account for variants within the same codon. I don't know weather it was a bug in the vcf file
output by GATKs unified genotyper or a bug in snpEff itself. But ultimately this code was used over snpEff because of this reason. 

I would reccomend using snpEff or a similar published program over this though. While it is accurate, it is also slow and fairly specific to the project I was working on in terms of the output that I wanted.
