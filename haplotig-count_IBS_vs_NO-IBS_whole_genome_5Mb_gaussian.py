import pandas as pd
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.backends.backend_pdf import PdfPages

# common imput sample name for multiple files
sample_name = (sys.argv[1])
path_file = (sys.argv[2])
# replace_sample_name = sample_name.replace(sample_name, '')

with PdfPages('{0}_{1}_whole_genome_gaussian.pdf'.format((path_file), (sample_name))) as pdf:
	# windows to capture each file
	windows = 500
	while windows <= 500:
	## create DataFrames
		count_df = pd.read_csv('{0}/whole_genome_haplo-tig_count_{1}_{2}kb_wind_gaussian_IBS.txt'.format((path_file), (sample_name), (windows)), delimiter='\t')
		count_df_no_ibd = pd.read_csv('{0}/whole_genome_haplo-tig_count_{1}_{2}kb_wind_gaussian_NON-IBS.txt'.format((path_file), (sample_name), (windows)), delimiter='\t')
		count_df_all = pd.read_csv('{0}/whole_genome_haplo-tig_count_{1}_{2}kb_wind_gaussian.txt'.format((path_file), (sample_name), (windows)), delimiter='\t')

		######### Set the number of plots
		fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(25,13), tight_layout=False)
		ax0,  ax1, ax2 = axes.flatten()

		## define variables for: ax0 = upper plot with raw data
		a = count_df['haplotig_count2'].astype(int)
		b = count_df_no_ibd ['haplotig_count2'].astype(int)
		c = count_df_all ['haplotig_count2'].astype(int)

		## make histograms for ax0
		ax0.hist(a, log=False, alpha=0.5, bins=range(min(a), max(a)), color='darkgreen', label='IBS')
		ax0.legend(prop={'size': 15})
		ax0.hist(b, log=False, alpha=0.5, bins=range(min(b), max(b)), color='darkblue', label='NO-IBS')
		ax0.legend(prop={'size': 15})
		ax0.grid(True)
		ax0.set_title('Haplo-tig count frequency IBS vs NO-IBS_gaussian [{0}:whole_genome {1}kb_wind]'.format((sample_name), (windows)), fontweight="bold", fontsize=15)
		ax0.set_ylabel('Haplo-tig count frequency', fontweight="bold", fontsize=15)
		ax0.set_xlabel('Haplo-tig count', fontweight="bold", fontsize=15)

		## make histograms for ax1
		ax1.hist(a, log=False, color='darkgreen', alpha=0.5, bins=range(min(a), max(a)), label='IBS')
		ax1.legend(prop={'size': 15})
		ax1.hist(b, log=False, color='darkblue', alpha=0.5, bins=range(min(b), max(b)), label='NO-IBS')
		ax1.legend(prop={'size': 15})
		ax1.grid(True)
		ax1.set_xscale('log')
		ax1.set_ylabel('Haplo-tig count frequency', fontweight="bold", fontsize=15)
		ax1.set_xlabel('Haplo-tig count [log]', fontweight="bold", fontsize=15)

		######### all_counts
		ax2.hist(c, log=False, color='darkred', alpha=1, bins=range(min(c), max(c)), label='all_regions')
		ax2.legend(prop={'size': 15})
		ax2.grid(True)
		ax2.set_xscale('log')
		ax2.set_ylabel('Haplo-tig count frequency', fontweight="bold", fontsize=15)
		ax2.set_xlabel('Haplo-tig count [log]', fontweight="bold", fontsize=15)

		# increase windows size
		windows += 100

		pdf.savefig()
