import os
import sys
import re 
import numpy as np 
import pandas as pd 
import subprocess

file_index = (sys.argv[1]) 				# chr_index.txt
file_length = (sys.argv[2]) 			# chr_sizes_jagger.txt
file_database = (sys.argv[3])			# jagger-cadenza_10x_concatenated.db 
chrom = (sys.argv[4])					# chromosome name
windows_resolution = int((sys.argv[5])) # windows resolution in bp

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


chrom_length = []
with open(file_length, 'r') as f_size:
	for line in f_size:
		chr_size = line.split("\t")
		# print(chr_size)
		if chrom == chr_size[0]:
			# print(chr_size[1])
			chrom_length = int(chr_size[1]) + windows_resolution
# print(chrom_length)

## Part 1: create template: this part takes less than a minute to run
# chrom_length = 500000 + 500000 # for testing
rep1M = 1 
window1M = 1 
rep5M = 1 
window5M = 1
start = 0
end = windows_resolution

stats_tenplate = open(file_database_name+'_'+chrom+'_haplo-tig_count_stats_score_'+str(windows_resolution)+'bp_windows.txt', 'w')
header = ('chr	rep	windows_1M	rep2	windows_5M	start	end' +'\n')
stats_tenplate.write(header)
stat_data = ''
while start < chrom_length and end <= chrom_length:
	stat_data = str(stat_data) + chrom +'\t'+ str(rep1M) +'\t'+ str(window1M) +'\t'+ str(rep5M) +'\t'+ str(window5M) +'\t'+ str(start) +'\t'+ str(end)+'\n'
	start = start + windows_resolution
	end = end +windows_resolution
	
	rep1M = rep1M + 1
	if rep1M == 11:
		window1M = window1M + 1
		rep1M = 1

	rep5M = rep5M + 1
	if rep5M == 51:
		window5M = window5M + 1
		rep5M = 1	
stats_tenplate.write(stat_data)
stats_tenplate.close()


# # Part 2. This part takes 2 minutes to run
# database_file = open(file_database)
# haplotig_size_file = open(file_database_name+'_'+chrom+'_haplo-tig_size-positions.txt', 'w')
# header_size = 'start_db	end_db	size'
# haplotig_size_file.write(header_size + ('\n'))
# for line in database_file:
# 	columns = line.strip('\n').split('\t')
# 	if int(columns[3]) >=chrom_indx and int(columns[3]) <(chrom_indx+10000000000): # check this part in other scripts since <= migh have been taking the 2 million posiitons as well
# 		# k_tig = columns[0]
# 		size = columns[1]
# 		# chromosome = columns[2]
# 		start = int(columns[3])-chrom_indx
# 		end = int(start) + int(size)
# 		output_data = str(start) + '\t' + str(end) + '\t' + str(size) + '\n'
# 		haplotig_size_file.write(output_data)
# haplotig_size_file.close()

# This part takes too long
nev_stat_data = pd.read_csv(file_database_name+'_'+chrom+'_haplo-tig_count_stats_score_'+str(windows_resolution)+'bp_windows.txt', delimiter='\t')
db_data = pd.read_csv(file_database_name+'_'+chrom+'_haplo-tig_size-positions.txt', delimiter='\t')

windows_1 = 0 
windows_2 = windows_resolution 
count_list = []
size_sum_list = [] 

while windows_2 <= chrom_length: 

	new_data_frame = db_data.loc[db_data['start_db'] >= windows_1, ['start_db']] < windows_2 
	size_sum = db_data.loc[(db_data['start_db'] >= windows_1) & (db_data['start_db'] < windows_2),  'size'].sum()

	windows_1 = windows_1 + windows_resolution 
	windows_2 = windows_2 + windows_resolution
	haplotig_count = new_data_frame.sum()
	count_list.append(haplotig_count[-1]) 
	size_sum_list.append(size_sum)

	nev_stat_data['haplotig_count']=pd.Series(count_list) 
	nev_stat_data['size_sum']=pd.Series(size_sum_list)
	
nev_stat_data['haplotig_score'] = nev_stat_data.haplotig_count * nev_stat_data.size_sum
# nev_stat_data.to_csv(file_database_name+'_'+chrom+'_haplo-tig_count_stats_score'+str(windows_resolution)+'bp_wind_count_only.txt', index=False, sep='\t')
nev_stat_data.to_csv(file_database_name+'_'+chrom+'_haplo-tig_count_stats_score_'+str(windows_resolution)+'bp_windows.txt', index=False, sep='\t')

# print("--- %s seconds ---" % round((time.time() - start_time)))



