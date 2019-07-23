# Haplotyope_report

This script reports haplotypes spanning 2 or more known SNPs per read. Simple counts are reported. Strongly inspired and influenced by LDx script by Alison Feder. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048588

# Usage:
LDx_haplo.py [-h] -i BAM -v VCF [-o OUT] [-d DEPTH] [-f MAF]
#arguments -i/--bam -v/--vcf are required

# Intended Application:
Viral quasispecies investigation - short range haplotype information is extremely useful for epitope quasispecies analysis. 

# Possible applications:
If you want to clarify whether any SNPs close enough to fall on the same read present multiple haplotypes in sequenced population.

# Note of Caution:
This script has not been optimized for speed. HBV genome is small (3Kb).
SNP files will need to be filtered and reformatted prior to usage. This script does not include VCF filtering.  
Filtering of low coverage haplotypes will be necessary.

# Future Improvements:
Allow user specification of region of interest. Currently, this can be done by pre-filtering SNPs of interest and reads that only map to the region of interest. 
At the moment this script is not equipped to handle data from more than 1 chromosome.
