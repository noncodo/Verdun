# Description of scripts 


Script | Description
------------ | -------------
main_commands.sh | Specific commands for basecalling, demultiplexing, consensus generation, variant filtering, etc. 
fix_medaka_fastas.sh | Shell script to process .vcf files and consensus genomes with ambiguous bases
rename_commands.sh | Shell script to help rename fastqs files output from guppy basecalling to sample-specific identifiers
vcf2fa.py | Python script to substitute ambiguous bases in consensus genome with most abundant variant when frequency is below ~90%
