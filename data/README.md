# Description of data files

File | Description
------------ | -------------
MN908947.ref.fasta.gz | SARS-CoV-2 reference genome (Wuhan-1)
all_medaka_variants_gt20_trim.vcf.gz | List of all filtered variants (supported by >20 reads)
medaka.fa.gz | Medaka-generated consensus genomes (c.f. ARTIC bioinformatics pipeline)
medaka_fixed.fa.gz | Medaka-generated and corrected with less abundant variant frequencies (see ../scripts/fix_medaka_fastas.sh) 
medaka_fixed_gt80pc.fa.gz | Subset of reads used for phylogenetic inference and lineage assignment 
nanopolish.fa.gz | Nanopore-generated consensus genomes (c.f. ARTIC bioinformatics pipeline)
negCtrls.fa.gz | Medaka-generated consensus genomes for negative controls (largely NNNNs)
pangolin.tgz.gz | Pangolineages for the above-mentioned fastas
phate_processed.tar | One-hot encoded clinical and variants data used to compute PHATE embeddings

