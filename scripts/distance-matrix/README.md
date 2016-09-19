# Distance matrices computed for following species:

+ k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pseudomonadales|f__Pseudomonadaceae|g__Pseudomonas|s__Pseudomonas_fluorescens 57 genomes
+ k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas|s__Stenotrophomonas_maltophilia 60 genomes
+ k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Burkholderiaceae|g__Burkholderia|s__Burkholderia_stagnalis 63 genomes
+ k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_pneumoniae 6659 genomes

# Algorithm to compute distance matrices:

+ Run ```subtrees_compute.py``` to generate directory of all genomes (in FASTA format) per species
+ Run PhyloPhlAn using genomes input directory
+ Convert full (concatenated) FASTA alignment to PHYLIP format using ```fasta2phylip.py```
+ Run Protdist on PHYLIP alignment
