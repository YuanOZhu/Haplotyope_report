#!/usr/bin/python
##########################################
#LDhaplo writer Yuan O Zhu: zhuy@gis.a-star.edu.sg
#Copyright (c) 2019, Yuan O Zhu, All rights reserved.
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#############################################
import os
import sys
import argparse
import numpy as np
#import pysam (not supported on aquila python2)
import subprocess #python2 uses .call not .run
import re
import time

##########################################
#Set up command line parameters and usage
##########################################
#usage usage: LDx_haplo.py [-h] -i BAM -v VCF [-o OUT] [-d DEPTH] [-f MAF]
#arguments -i/--bam -v/--vcf are required
#eg. python LDx_haplo.py -i sample.sorted.grouped.realign.recal.bam -v sample.filtered.vcf -o out.txt
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--bam', help='sorted bam files', required=True)
parser.add_argument('-v', '--vcf', help='filtered vcf with positions in 2nd column', required=True)
parser.add_argument("-o", "--output", type=argparse.FileType('w'), help="Directs the output to a file of your choice", default="haplo.out")
parser.add_argument('-c','--counts', help='minimum reads a haplotype must be observed to be reported', default=2)
parser.add_argument('-f','--maf', help='only report haplotypes above this frequency', default=0.05)
parser.add_argument('-b','--bases', help='minimum read length', default=100)
parser.add_argument('-p','--palign', help='minimum % of read aligned per read', default=0.7)
parser.add_argument('-q','--mapq', help='minimum read map quality', default=30)
args = parser.parse_args()
output_file=args.output
output_log="LDx.log"
log=open(output_log, "a+")

############################################################################
#Initialize some variables to keep track of things we'll want at the end
############################################################################
#List of known SNP locations [pos1,pos2,pos3]
SNPs=[]
#Dictionary of haplotype counts
#key=pos1,pos2,pos3-base1,base2,base3 value=counts observed
haplo_count={}
#corrected dictionary
haplo_new={}
#tracks total coverage over every haplotype position
haplo_sum={}
#for new haplotypes generated during iter
#haplo_count_2={}

#####################################################################################################################
#SNP processing
#To allow for flixibility in vcf format, vcf files are expected to be pre-filtered for quality and allele frequency
#Input vcf files need only have a second column stating the position of identified SNPs
#####################################################################################################################
SNPs = np.genfromtxt(args.vcf, delimiter="\t", dtype=str, usecols=1)
#print(SNPs)


#####################################################################################################################
#define subroutines required by main execution code
#identifies reads covering >2 SNP positions and stores SNP position+haplotype in dictionary
#####################################################################################################################
#passMap takes in cigar string and returns 1 if read passes length and % mapped filters
def passMap (mapcigars):
    booleanreturn = 0
    mapping = re.split('(\d+)',mapcigars) #breakdown mapping cigar string
    numbers = map(int,mapping[1::2]) #list of cigar numbers
    texts = mapping[2::2] #corresponding list of cigar alphabets
    mappedlength = 0
    stringlength = sum(numbers) #print ("stringlength "+ str(stringlength))
    if (len(texts) < 4): #filters out reads with too many indels
        for i in range(len(texts)):
            if (texts[i] == "M"):
                mappedlength += numbers[i] #print ("mappedlength "+ str(mappedlength))
        if (stringlength == 0):
            booleanreturn = 0
        else:
            percentmapped=float(mappedlength)/float(stringlength) #print("percentmapped "+str(percentmapped))
        if ((stringlength > args.bases) and (percentmapped > args.palign)):
            booleanreturn = 1
        else:
            booleanreturn = 0
    else:
        booleanreturn = 0        
    return booleanreturn

#extractHaplo takes nucleotide string and mapped coordinates on reference
#returns base calls on read that map to known SNPs if no. SNPs>1 
def extractHaplo(ref_pos, ref_end, segment_string):
    haplo_posi = "" #string to store all SNP positions covered
    haplo_sequ = "" #string to store all bases called at SNP positions covered
    snp_count = 0 #tracks number of SNPs covered
    for x in range(int(ref_pos), (int(ref_end))): #for every mapped position on reference
        #print ("ref_position "+str(x))
        str_pos = int(x) - int(ref_pos) #find corresponding position on string
        #print ("string_position "+str(str_pos))
        str_base = segment_string[str_pos] #find corresponding base call on string
        #print ("string_base "+str_base)
        if str(x) in SNPs: #if mapped position is a known SNP position
            #print("Known SNP:"+str(x)+" string_position:"+str(str_pos)+" string_base:"+str(str_base))
            if (snp_count == 0): #attach position and base information to strings
                haplo_posi = str(x)
                haplo_sequ = str(str_base)
            else:
                haplo_posi = haplo_posi + "," + str(x)
                haplo_sequ = haplo_sequ + "," + str(str_base)
            snp_count += 1
    if (snp_count > 1): #if more than 1 SNP in read
        #define key value of positions+genotype
        key=haplo_posi+"|"+haplo_sequ
        #add count to haplo_sum
        if haplo_posi in haplo_sum:
            haplo_sum[haplo_posi]+=1
        else:
            haplo_sum[haplo_posi]=1
        #add count to haplo_count
        if key in haplo_count:
            haplo_count[key]+=1
        else:
            haplo_count[key]=1

#haploCut determines whether read covers >2 SNPs and stores haplotype if so
def haploCut (startpos, cigar, seq):
    read_pos = 0 #tracks start of read segment under scan base-1
    ref_pos = startpos #tracks start position on reference exact base
    ref_end = startpos #tracks end of segment on reference exact base
    haplo_seq = "" #string to store haplo_sequence for read
    haplo_pos = "" #string to store haplo_positions for read
    mapping = re.split('(\d+)',cigar) #breakdown mapping cigar string    
    numbers = map(int,mapping[1::2]) #list of cigar numbers
    texts = mapping[2::2] #list of cigar alphabets
    discards = ["N","X","="] #list of cigar alphabets to discard read on (for simplicity)
    if not (any(i in discards for i in texts)): #ignore reads with N/X/= in cigar
        for i in range(len(texts)):
            if (texts[i] == "M"): #continue with matched bases
                ref_end = int(ref_pos) + int(numbers[i]) -1
                #print("Matched read pos"+str(read_pos)+" "+str(i)+" "+str(numbers[i]))
                #print("Matched ref range"+str(ref_pos)+" "+str(ref_end))
                segment_string = seq[read_pos:(read_pos+numbers[i]+1)]
                #print("Matched segment "+segment_string)
                extractHaplo(ref_pos, ref_end, segment_string)
                ref_end += 1
                ref_pos = ref_end
                read_pos = int(read_pos) + int(numbers[i])
            elif (texts[i] == "S"): #skip positions on read for soft trim
                read_pos = int(read_pos) + int(numbers[i]) #print("Soft trim| read_pos:"+str(read_pos))
            elif (texts[i] == "D"): #insert positions on read for deletion
                ref_end = int(ref_pos) + int(numbers[i]) -1 #print("Deletion| ref_pos:"+str(ref_pos)+" ref_end:"+str(ref_end))
            elif (texts[i] == "I"): #skip positions on read for insertion
                read_pos = int(read_pos) + int(numbers[i]) #print("Insertion| read_pos:"+str(read_pos))
            #elif ((texts[i] == "H") or (texts[i] == "P")): #ignore if hard trim or padded
            #print ("Loop ref_pos:"+str(ref_pos)+" ref_end:"+str(ref_end))

def updateHaplo():
    for y in haplo_count:
        #print(y+"\t"+str(haplo_count[y]))
        yout = y.split("|")
        haplos=yout[0].split(",")
        seq=yout[1].split(",")
        length=len(seq)
        if(length > 2):
            for i in range(2,length):#for each winsize from 2 to length of haplo-1
                winsize=i
                for j in range(0,length-winsize+1):#for each position
                    temp_seq = seq[j:(j+winsize)]
                    temp_pos = haplos[j:(j+winsize)]
                    input_seq = ",".join(temp_seq)
                    input_pos = ",".join(temp_pos)
                    temp_id = input_pos + "|" + input_seq
                    if input_pos in haplo_sum:
                        haplo_sum[input_pos] += haplo_count[y]
                    else:
                        haplo_sum[input_pos] = haplo_count[y]
                    if(temp_id in haplo_new):
                        haplo_new[temp_id] += haplo_count[y]
                    else:
                        haplo_new[temp_id] = haplo_count[y]

                        
#print haplo counts if haplotype is observed > user defined times (default=2)
def finalOut(dictN,dictSum):
    for x in dictN:
        xout = x.split("|")
        if (int(dictN[x]) > int(args.counts)):
            maf=float(dictN[x])/float(dictSum[xout[0]])
            if(float(maf) > float(args.maf)):
                if(float(maf) < (1-float(args.maf))):
                    output_file.write(xout[0]+"\t"+str(dictSum[xout[0]])+"\t"+xout[1]+"\t"+str(dictN[x])+"\t"+str(maf)+"\n")
                
                        
############################## THIS IS THE MAIN EXECUTION CODE ####################################           
#identifies reads covering >2 SNP positions and stores SNP position+haplotype in dictionary 
###################################################################################################

def main ():
    #convert bam to sam for parsing
    bamfile=args.bam
    samfile=bamfile+".sam"
    #subprocess.check_call("samtools view -sB -q args.mapq "+bamfile+" > "+samfile, shell=True)
    subprocess.check_call("samtools view -sB "+bamfile+" > "+samfile, shell=True)
    
    #read in sam file, identify pass quality reads and identify snp haplotypes
    with open(samfile) as fp:
        line = fp.readline()
        while line:
            stripped = line.rstrip('\n')
            data = line.split("\t")
            if((data[5] != "*") and (int(data[4]) >= args.mapq)): #for mapped reads with mapq>=30
                #print("\n\nMAIN "+stripped)
                startbase=data[3] #first base match position
                mapcigars=data[5]
                sequence=data[9] #DNA sequence
                if (passMap(mapcigars) == 1): #matched bases must be continuous and > 60% of read #print ("passed "+line)
                    haploCut(startbase, mapcigars, sequence)
            line = fp.readline()
    fp.close()

    updateHaplo() #for every long haplo of n>2, add counts for all encompassed short haplos length 2 to n-1

    #order updated haplo dictionary
    ordered_haplo=sorted(haplo_count)

    #print file header
    output_file.write("Haplo_pos\tHaplo_coverage\tHaplo_genotype\tGeno_count\tGeno_freq\n")

    #merge overlapping haplotypes
    final_haplos = {}
    for key in ordered_haplo:
        if key in haplo_new:
            final_haplos[key]=haplo_count[key]+haplo_new[key]
            haplo_count.pop(key)
            haplo_new.pop(key)
        
    finalOut(final_haplos,haplo_sum)
    finalOut(haplo_new,haplo_sum)
    finalOut(haplo_count,haplo_sum)

    
###################################################################################################
############################################ RUN CODE #############################################
###################################################################################################

start_time = time.time()
main()
log.write("\nLDx.py run on file " + args.bam + " completed in --- %s seconds ---" % (time.time() - start_time))
