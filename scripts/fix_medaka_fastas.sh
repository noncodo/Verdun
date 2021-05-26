#!/bin/bash
################################################################################################
## merge and parse vcfs into variant files to include variant allele frequencies below 90%
## (c) Martin A. Smith 2021, University of Montreal
##
##  Usage : ./fix_medaka_fastas.sh artic_output_directory
##
################################################################################################

INDIR=$1

#gunzip ${INDIR}/*vcf.gz ## uncomment this for first time
echo -ne "[*] Parsing variant files..."
for file in ${INDIR}/*.pass.vcf ${INDIR}/*.fail.vcf ; do 
	STAT=${file%*.vcf}
	STAT=${STAT#*.}
	ID=${file%%.*} 
	ID=${ID#*/}

   	# Check if file already exists, meaning the first (pass/fail) vcf was processed
	# This will help overwrite past files if re-running to avoit appending >2 times
	if [[ ! -e medaka2/${ID}.filtVars.tsv ]]; then
   		FIRST=1
	else 
   		FIRST=0
	fi
	
	# clean up .vcfs into easy to parse .tsv format
	if [[ ${FIRST} -eq 1 ]]; then 
		grep -v ^# ${file} | sed -e 's/;/\t/g' -e 's/DP=//g' -e 's/AC=//g' -e 's/,/\t/g' -e 's/AM=//g' -e 's/AQ=//g' |\
			awk -v id=$ID -v st=$STAT 'OFS="\t"{print id,$2,$3,$4,$5,$6,st,$8,$9,$10,$10/($9+$10),$11,$15}' > ${INDIR}/${ID}.filtVars.tsv
		FIRST=0
	else
		grep -v ^# ${file} | sed -e 's/;/\t/g' -e 's/DP=//g' -e 's/AC=//g' -e 's/,/\t/g' -e 's/AM=//g' -e 's/AQ=//g' |\
			awk -v id=$ID -v st=$STAT 'OFS="\t"{print id,$2,$3,$4,$5,$6,st,$8,$9,$10,$10/($9+$10),$11,$15}' >> ${INDIR}/${ID}.filtVars.tsv
	fi	
done 
echo " DONE"

################################################################################################
# now let's use custom filters to edit the reference genome based on selected variants 
# By default, the medaka artic pipeline treats variants with <90% VAF as 'N'
# we want to lower that threshold tp 50% to pick up the most abundant variant in the consensus
echo "[*] Fixing consensus fasta files"
echo
for file in ${INDIR}/*.filtVars.tsv ; do
	# tested on python 3.7+
	echo -e '\e[1A\e[K     '$file
	python vcf2fa.py -f MN908947.ref.fasta -v ${file} -c ${file%*.filtVars.tsv}.coverage_mask.txt
done
echo "[!] All done :-)"
