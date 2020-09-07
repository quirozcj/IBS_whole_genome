import pandas as pd
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter, FixedLocator
from contextlib import ExitStack
import matplotlib
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)

path_file = (sys.argv[1])
sample_name = (sys.argv[2])
components_no = int((sys.argv[3]))
chromosome_id = (sys.argv[4])
WINDOWS = (sys.argv[5])
FILTER = (sys.argv[6])
chr_length = (sys.argv[7])

PR_haplotig_position = pd.read_csv('{0}precision_recall_whole_genome_haplo-tig_count_{1}_{2}_kb_wind_filter_{3}_gauss_comp{4}.txt'.format((path_file), (sample_name), (WINDOWS), (FILTER),(components_no)), delimiter='\t')

plt.figure(figsize=(30, 8))
fig, ax = plt.subplots(figsize=(30, 8))
plt.suptitle('Precision-Recall: '+chromosome_id+'_'+sample_name+'_'+str(WINDOWS)+'_windows_'+str(FILTER)+'_filter', fontweight="bold", fontsize=20)
plt.grid(True, color='gray', linestyle='--', linewidth=0.5)


FP_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chromosome_id) & (PR_haplotig_position['FP'] == 1)]
x = FP_chrom_frames['start']
y = FP_chrom_frames['haplotig_count']
plt.scatter(x,y, c = "blue", alpha = 1, s = 10,label='FP')
plt.legend(prop={'size': 10})


FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chromosome_id) & (PR_haplotig_position['FN'] == 1)]
a = FN_chrom_frames['start']
b = FN_chrom_frames['haplotig_count']
plt.scatter(a,b, c = "darkred", alpha = 1, s = 10, marker='>',label='FN')
plt.legend(prop={'size': 10})

FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chromosome_id) & (PR_haplotig_position['TP'] == 1)]
c = FN_chrom_frames['start']
d = FN_chrom_frames['haplotig_count']
plt.scatter(c,d, c = "darkgreen", alpha = 1, s = 10, marker='<',label='TP')
plt.legend(prop={'size': 10})

# FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chromosome_id) & (PR_haplotig_position['gauss_mx_IBS'] == 1)]
# g = FN_chrom_frames['start']
# h = FN_chrom_frames['haplotig_count']
# plt.scatter(g,h, c = "turquoise", alpha = 1, s = 10,,label='IBS')
# plt.legend(prop={'size': 10})

FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chromosome_id) & (PR_haplotig_position['TN'] == 1)]
e = FN_chrom_frames['start']
f = FN_chrom_frames['haplotig_count']
plt.scatter(e,f, c = "orange", alpha = 1, s = 2, marker='*',label='TN')
plt.legend(prop={'size': 10},markerscale=2,borderaxespad=0.5,handletextpad=.2) # label legen position

######## change this part for a DataFrame pandas#######

chr_length_list = []
with open(chr_length, 'r') as f_size:
	for line in f_size:
		chromosome_size = line.split("\t")
		if chromosome_id == chromosome_size[0]:
			chr_length_list = int(chromosome_size[1])

chr_length_list=chr_length_list+5000000
chr_length_list_labels=int(chr_length_list/1000000)
plt.xlim(0,chr_length_list)
x_labels=list(range(0,chr_length_list,10000000))
x_ticklabels=list(range(0,chr_length_list_labels,10))
plt.xticks(x_labels, x_ticklabels, rotation=45 )

######## change this part for a DataFrame pandas#######

# plt.xlim(0,750000000)
# x_labels=list(range(0, 750000000, 10000000)) 
# x_ticklabels=list(range(0,750,10)) 
# plt.xticks(x_labels, x_ticklabels, rotation=45 )
plt.yscale('symlog', base=10)
# plt.yticks([0,1,2,5,10,30,50,100,200])
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.set_ylabel('Haplo-tig count [Log]', fontweight = "bold", fontsize = 15)
ax.set_xlabel('Chromosome position [Mbp]', fontweight = "bold", fontsize = 15)

# plt.show()
plt.savefig(chromosome_id+'_'+sample_name+'_'+str(WINDOWS)+'_windows_'+str(FILTER)+'_filter''.pdf',dpi = 75, transparent = False, bbox_inches='tight')
	# plt.close()