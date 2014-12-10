#!/usr/bin/env python

"""
Parse Exon and Mask non-exonic regions in the repeat and polymorphism masked sequences

Description:    For each of the repeat and polymorphism masked sequences, the correpsonding exon coordinates are extracted from ZmB73_5a.59_WGS.gff.gz.tmp file.
                This is based on their chr no and start and stop positions. Regions outside of the exon coordinates are masked out as X.
                This non-exon masked sequence can be used for designing LR PCR primers.

"""


import time                     # timer function
import glob                     # function to read all file names in current directory
import datetime
now = datetime.datetime.now()
date_month_year = str(now.month)+'_'+ str(now.day)+ '_' +str(now.year)


# Time to run the code: start timer
t0 = time.time()







# Read all file names in sequences directory, that ends with.txt
file_names =  (glob.glob('sequences/*.txt') )
print file_names
for name in file_names:
                print name[10:]
                cols=name.split('_')
                chr_no = int(cols [1])
                start = int(cols [2])
                stop = int(cols [2]) + int(cols [3])-1

 

# Input file (ZMDB working set annotations)
                with open (name) as sequence_data:
                                sequence = sequence_data.read()
                                seq_end_coord = len(sequence)
                
                exon_start_list = []
                exon_stop_list = []               
                #with open ('test_ZmB73_5a.59_WGS.gff.gz.tmp') as input_data:
                with open ('ZmB73_5a.59_WGS.gff.gz.tmp') as input_data:
                                lines = input_data.readlines()
                                for line in lines:
                                                cols1=line.split('\t')
                                                if cols1[0].isdigit() and ((cols1[2]) == 'exon') and (int(cols1[0]) == chr_no) and (int(cols1[3]) >= start) and (int(cols1[4]) <= stop):                                                            
                                                #if cols1[0].isdigit() and ((cols1[2]) == 'exon') and (int(cols1[0]) == chr_no) and (int(cols1[3]) >= start) and (int(cols1[4]) <= stop):
                                                                exon_start = int(cols1[3]) - start
                                                                exon_stop = int(cols1[4]) - start
                                                                exon_start_list.extend([exon_start])
                                                                exon_stop_list.extend([exon_stop])


# Mask non exonic regions in the sequence               
                with open (name) as sequence_data:
                                sequence = sequence_data.read()
                                masked_sequence = "X" * len(sequence)
                                for i, j in zip(exon_start_list, exon_stop_list):
                                                #print sequence[i-1:j]
                                                if i > 0:
                                                                masked_sequence =  (masked_sequence[0:i-1]+sequence[i-1:j]+masked_sequence[j-1:-1] )
                                                if i == 0:
                                                                masked_sequence =  (sequence[i:j]+masked_sequence[j-1:-1] )
                                                                

                out_put = masked_sequence
                print len(sequence)
                print len(masked_sequence)
                if len(masked_sequence) != len(sequence):
                                out_put = "Error! Check your annotation file"
                print (out_put)
                
# Write to output.txt file

                with open('Non_exon_masked'+ date_month_year + '_' + name[10:], 'w') as output_data:
                                output_data.write(out_put)


# Time to run the code: end timer
t1 = time.time()
total = t1-t0
print (total)

