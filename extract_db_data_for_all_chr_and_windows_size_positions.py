import os
import sys
import re 
import numpy as np 
import pandas as pd 
import subprocess

file_index = (sys.argv[1]) 				# chr_index.txt
file_database = (sys.argv[2])			# jagger-cadenza_10x_concatenated.db 
chrom = (sys.argv[3])					# chromosome name

# define name for output files
file_database_name = file_database.replace('_concatenated.db', '')

# index_file
chrom_indx = []
with open(file_index, 'r') as f_index:
	for line in f_index:
		chr_indx = line.split("\t")
		if chrom == chr_indx[0]:
			chrom_indx = int(chr_indx[1])
# print(chrom_indx)

# ## Part 2. This part takes 2 minutes to run
database_file = open(file_database)
# if not os.path.exists(file_database_name+'_'+chrom+'_haplo-tig_size-positions.txt'):
haplotig_size_file = open(file_database_name+'_'+chrom+'_haplo-tig_size-positions.txt', 'w')
header_size = 'start_db	end_db	size'
haplotig_size_file.write(header_size + ('\n'))
for line in database_file:
	columns = line.strip('\n').split('\t')
	if int(columns[3]) >=chrom_indx and int(columns[3]) <(chrom_indx+10000000000): # check this part in other scripts since <= migh have been taking the 2 million posiitons as well
		# k_tig = columns[0]
		size = columns[1]
		# chromosome = columns[2]
		start = int(columns[3])-chrom_indx
		end = int(start) + int(size)
		output_data = str(start) + '\t' + str(end) + '\t' + str(size) + '\n'
		haplotig_size_file.write(output_data)
haplotig_size_file.close()



