#!/usr/bin/python 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", required=True,
                       help="Fasta reference file")
parser.add_argument("-v", "--var", required=True,
                       help="tab-delineated file of variants")
parser.add_argument("-c", "--mask", required=True, 
                       help="Coverage mask .txt file to print out number of ambiguopus bases")
parser.add_argument("-t", "--vaf", type=float, default=0.5,
                       help="variant allele threshold for inclusion into fasta (def. 0.5)")
parser.add_argument("-d", "--depth", type=int, default=20,
                       help="minimum read depth to include a variant (def. 20)")
args = parser.parse_args()

vaf = args.vaf
dic = dict()

# import variants as dictionary
with open(args.var, 'r') as v:
 for line in v:
 ###### Output of fix_medaka_fastas.sh. LOOKS LIKE: 
 # BA8 18361 . CA  TG  32.37 fail  2 1 1 0.5 0 37.75
 # BA8 21408 . TC  T 78.65 fail  62  32  17  0.346939  13  10.71
 # BA8 199 . G GT  500.0 pass  2491  350 1910  0.845133  231 12.61
 # BA8 241 . C T 500.0 pass  2490  9 2280  0.996068  201 15.41
 # BA8 3037  . C T 500.0 pass  665 3 600 0.995025  62  12.21
  if line[0] != '#':
    s=line.strip("\n").split()
    dic[s[1]] = s[3] + " " + s[4] + " " + s[10] + " " + s[7]
    #print(dic)

# get the consensus genome
header=""
seq=""
file=args.var.split("/")[-1]
name=file.split(".")
with open(args.fasta, 'r') as f:
 for line in f: 
  if line[0] == ">":
   header=">"+name[0]+" medaka_fix\n"
  else:
   seq=seq+line.strip()

# fix the medaka consensus genome
# this only works because dict. outputs a sorted index
last_pos=0
new_seq=""
new_seq2=""
for v in dic:
  var=dic[v].split() 
  # apply custom thresholds (min VAF and coverage depth per variant)
  if float(var[2]) >= vaf and int(var[3]) > args.depth:
    var_size=len(var[0])-len(var[1]) 
    if var_size > 0 : # DELETION
      #print("Deletion", var_size, var[0], var[1], var[2] )
      new_seq2=new_seq+seq[last_pos:int(v)-1] + var[1]
      last_pos=int(v)+var_size # makes sure the deleted chars are gone
    elif var_size < 0 : # INSERTION
      #print("Insertion ", var_size , var[0], var[1], var[2] )
      new_seq2=new_seq+seq[last_pos:int(v)-1] + var[1]
      last_pos=int(v)
    else : 
      new_seq2=new_seq+seq[last_pos:int(v)-1] + var[1]
      last_pos=int(v)
    new_seq=new_seq2  
new_seq2=new_seq+seq[last_pos:]
new_seq=new_seq2  

#write the fixed fasta
outfile=args.var.rsplit('.',2)[0] + ".consensus_fixed.fasta"
counter=int()
with open(outfile, 'w') as o:
 o.write(header)
 for i in range(0,len(new_seq),60):
 	o.write(new_seq[i:i+60] + "\n")

