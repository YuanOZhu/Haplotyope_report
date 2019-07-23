# Haplotyope_report

This script reports haplotypes spanning 2 or more known SNPs per read. Simple counts are reported. Influenced by LDx script by Alison Feder. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048588

The output file will read like this
Haplo_pos         Haplo_coverage  Haplo_genotype  Geno_count	Geno_freq\
95,118            311185          A,C             171506      0.551138390347\
95,118            311185          A,T             115967      0.37266256407\
95,118            311185          G,T              19861      0.063823770426\
95,118,159        143135          A,C,A            35268      0.246396758305\
95,118,159        143135          A,C,G            45613      0.318671184546\
95,118,159        143135          A,T,A            49357      0.344828308939\
95,118,159        143135          G,T,A             7162      0.05003667866\
95,118,159,189     93117          A,C,A,T          22917      0.246109732917\
95,118,159,189     93117          A,C,G,T          29841      0.320467798576\
95,118,159,189     93117          A,T,A,T          32109      0.344824253359

Haplo_pos: positions of SNPs in haplotype\
Haplo_coverage: number of reads that covered this exact haplotype\
Haplo_genotype: genotype being reported\
Geno_count: how many reads covering the haplotype reported this genotype\
Geno_freq: Geno_count/Haplo_coverage

# Usage:
LDx_haplo.py [-h] -i BAM -v VCF [-o OUT] [-d DEPTH] [-f MAF]\
#arguments -i/--bam -v/--vcf are required

# Intended Application:
Viral quasispecies investigation - short range haplotype information is extremely useful for epitope quasispecies analysis.

# Possible applications:
If you want to clarify whether any SNPs close enough to fall on the same read present multiple haplotypes in sequenced population.

# Note of Caution:
This script has not been optimized for speed. HBV genome is small (3Kb).
SNP files will need to be filtered and reformatted prior to usage. This script does not include VCF filtering.

# Future Improvements:
Allow user specification of region of interest. Currently, this can be done by pre-filtering SNPs of interest and reads that only map to the region of interest.
At the moment this script is not equipped to handle data from more than 1 chromosome.
