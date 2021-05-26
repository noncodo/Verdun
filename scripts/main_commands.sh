#!/bin/bash


GUPPY=4.4.1

########################
#  Basecalling with guppy --version 4.4.1
find /data/nanopore/Verdun/verd? /data/nanopore/Verdun/verd?? -type d -name "fast5" | while read f5 ; do
  run=$( echo $f5 | cut -d '/' -f 5 ) 
  /home/apps/ont-guppy/ont-guppy-${GUPPY}/bin/guppy_basecaller guppy_basecaller \
	 	--config dna_r9.4.1_450bps_hac.cfg \
	 	--gpu_runners_per_device 272 \
	 	--chunk_size 2000 \
	 	--chunks_per_runner 1024 \
	 	-x "cuda:0 cuda:1 cuda:2 cuda:3" \
	 	--disable_pings \
	 	--input_path ${f5} --save_path /home/martin/VERDUN/${run}/fastq/guppy-${GUPPY}
done

########################
#  Demultiplexing with guppy --version 4.4.1
cd /home/martin/VERDUN/
for RUN in verd{1..14} ; do \
	echo "[*] demuxing "${RUN} 
	/home/apps/ont-guppy/guppy-${GUPPY}/bin/guppy_barcoder \
 		-x "cuda:0 cuda:1 cuda:2"   \
 		--require_barcodes_both_ends \
 		-i ${RUN}/fastq/ \
 		-s ${RUN}/fastq/demuxed/guppy-${GUPPY}/ \
 		--arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
done

#########################
# Renaming the fastqs
sed 's/VERD3_2/VERD3/g' sample_batches.txt |\
 awk 'NF>1' | tail -n +2 | sed 's/VERD/verd/g' |\
 awk 'OFS="\t"{ printf "cat "$2"/fastq/demuxed/guppy-4.4.1/";  
 	if (length($3)==1) printf "barcode0"$3"/fastq_runid_* " ; 
 	else printf "barcode"$3"/fastq_runid_*" ; print "> fastq_demux/"$1".fastq"}' > rename_commands.sh

#########################
# Collating the negative controls
mkdir negCtrl
awk '{ printf "cat "$1"/fastq/demuxed/guppy-4.4.1/barcode"$2"/fastq_runid_* > negCtrl/negCtrl_"NR".fastq\n"}' negCtrl.ids > negCtrl/collate_fastqs.sh
chmod 755 collate_fastqs.sh && ./collate_fastqs.sh 

#########################
# Size filtering
for file in ./fastq_demux/*fastq ; do  
  awk '{ if (NR%4==0) print $1 ; else printf $1"\t" }' $file | awk '{ len=length($2) ; if (len <= 700 && len >= 400 ) print $1,"\n"$2"\n"$3"\n"$4}' > fastq_demux_filt/${file##*/} 
done
for file in ./negCtrl/fastq/*fastq ; do  
  awk '{ if (NR%4==0) print $1 ; else printf $1"\t" }' $file | awk '{ len=length($2) ; if (len <= 700 && len >= 400 ) print $1,"\n"$2"\n"$3"\n"$4}' > ./negCtrl/fastq_filt/${file##*/} 
done

#########################
# Running the ARTIC pipeline
# Version 1.2.1 (conda)
# https://github.com/artic-network/fieldbioinformatics
###### NANOPOLISH VERSION 
conda activate artic
i=0
N=4
for file in fastq_demux_filt/*fastq ; do  
  ID=${file##*/}  
  ID=${ID%*.fastq}
  ((i=i%N)) 
  ((i++==0)) && wait 
  artic minion --normalise 2000 --threads 4 --scheme-directory /home/apps/artic-ncov2019/primer_schemes --read-file $file --fast5-directory /scratch/VERDUN/fast5/ --sequencing-summary  /home/martin/VERDUN/all_seq_sum.txt nCoV-2019/V3 $ID 2>/dev/null & 
done
###### MEDAKA VERSION 
#this activates CPU version of medaka because I messed up TensorFlow and the GPU version was a pain to get working
for file in ./fastq_demux_filt/*fastq ; do
  NAME=${file##*/} 
  NAME2=${NAME%*.fastq}
  artic minion --medaka --normalise 2000 --threads 24 --scheme-directory /home/apps/artic-ncov2019/primer_schemes --read-file $file nCoV-2019/V3 ${NAME2} 
done

#########################
# Fixing medaka consensus 
./fix_medaka_fastas.sh ./medaka/

#########################
# Running Pangolin
conda activate pangolin 
pangolin -v
# pangolin 2.1.10
pangolin -t 24 nanopolish_consensii.fa -o pangolin/nanopolish
pangolin -t 24 consensii.fa -o pangolin/medaka
pangolin -t 24 all_medaka_fixed.fasta -o pangolin/medakaFix
pangolin -t 24 negCtrls.fa -o ./pangolin/negCtrl
cd pangolin && tar -c m* n* | gzip > pangolin.tgz