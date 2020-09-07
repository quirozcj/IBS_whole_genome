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
# replace_sample_name = sample_name.replace(sample_name, '')

with PdfPages('{0}_FP_rcll_genome_pos.pdf'.format(sample_name)) as pdf:
	# for WINDOWS in range(100, 1100, 100):
	for WINDOWS in list([100, 200, 500, 1000]):	
		for FILTER in list([200, 600, 1000, 2000]):

			# Read the recall-position file
			PR_haplotig_position = pd.read_csv('{0}rcll_{1}_{2}_wind_{3}_flt_cp{4}.txt'.format((path_file), (sample_name), (WINDOWS), (FILTER),(components_no)), delimiter='\t')

			#################################################################
			## create figure:
			fig, axes = plt.subplots(nrows=7, ncols=3, sharey=False, sharex=True, figsize=(25,10))
			# fig.tight_layout()
			plt.subplots_adjust(top = 0.95, left=None, right = None, wspace=0.08, hspace=None)


			### Fugure main title:
			fig.suptitle('Precision-Recall FP vs FN: {0}_{1}_kb_windows_{2}_filter_components_{3}'.format((sample_name), (WINDOWS), (FILTER), (components_no)),fontweight="bold", fontsize=15)

			# chromosome limits
			plt.xlim(0,900000000)
			x_labels=list(range(0, 900000000, 50000000)) 
			x_ticklabels=list(range(0,900,50)) 
			plt.xticks(x_labels, x_ticklabels, rotation=45 )


			# # haplotig-count limits 
			# plt.ylim(0, 205)
			# y_labels=list(range(0, 205, 50))
			# y_ticklabels=list(range(0,205,50))
			# plt.yticks(y_labels, y_ticklabels)

			# Figure axes array: Note sublpots uses array for flatten method, it does not work with lists
			axes_list =[]
			for ax in range(1,22):
				axes_list.append('ax'+str(ax))
			axes_array = np.array(axes_list)
			axes_array = axes.flatten()

			# get labels title for y axes sub plots
			for y_axis in range (0, 19, 3):
				axes_array[y_axis].set_ylabel('haplotig-count', fontweight = "bold", fontsize = 8)

			# plot each chromosome in a loop for each subplot
			axis = 0
			for num in range(1,8):
				for chrom in list('ABD'):
					chrom_name = ('chr'+str(num)+(chrom))
			        # filter for chromosome name and False Positives olny
					FP_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chrom_name) & (PR_haplotig_position['FP'] == 1)]
					x = FP_chrom_frames['start']
					y = FP_chrom_frames['haplotig_count']

					axes_array[axis].set_xticklabels(x_ticklabels, rotation=45, fontsize = 10)
					axes_array[axis].scatter(x,y, c = "darkblue", alpha = 0.3, s = 3, label='FP')
					axes_array[axis].grid(True)
					axes_array[axis].set_xlabel(chrom_name +': positions [Mbp]', fontweight = "bold", fontsize = 8)
					axes_array[axis].set_yscale('symlog',base=10)
					# plt.legend(prop={'size': 5})
					

					FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chrom_name) & (PR_haplotig_position['FN'] == 1)]
					a = FN_chrom_frames['start']
					b = FN_chrom_frames['haplotig_count']

					axes_array[axis].set_xticklabels(x_ticklabels, rotation=45, fontsize = 10)
					axes_array[axis].scatter(a,b, c = "darkred", alpha = 0.3, s = 3, marker='>',label='FN')
					axes_array[axis].grid(True)
					axes_array[axis].set_xlabel(chrom_name +': positions [Mbp]', fontweight = "bold", fontsize = 8)
					axes_array[axis].set_yscale('symlog',base=10)
					axes_array[axis].set_yticks([0,1,3,10,40]) # specify the thicks we want only
					# plt.legend(prop={'size': 5})

					FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chrom_name) & (PR_haplotig_position['TP'] == 1)]
					c = FN_chrom_frames['start']
					d = FN_chrom_frames['haplotig_count']
					axes_array[axis].scatter(c,d, c = "darkgreen", alpha = 0.3, s = 3, marker='<',label='TP')
					# plt.legend(prop={'size': 5})


					FN_chrom_frames = PR_haplotig_position[(PR_haplotig_position['chr'] == chrom_name) & (PR_haplotig_position['TN'] == 1)]
					e = FN_chrom_frames['start']
					f = FN_chrom_frames['haplotig_count']
					axes_array[axis].scatter(e,f, c = "orange", alpha = 0.3, s = 3, marker='*',label='TN')
					plt.legend(prop={'size': 5},markerscale=2,borderaxespad=0.5,handletextpad=.2) # label legen position



					# axes_array[axis].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
					# axes_array[axis].yaxis.set_major_formatter(ScalarFormatter())
					# axes_array[axis].yaxis.set_minor_formatter(NullFormatter())
					# axes_array[axis].set_ylim(bottom = -0.5, top=None)
					# plt.ylim(bottom=0)
					# axes_array[axis].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
					# axes_array[axis].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
					# axes_array[axis].yaxis.set_minor_locator(ticker.MaxNLocator(5))
					# axes_array[axis].yaxis.set_major_locator(ticker.MaxNLocator(5))
					# axes_array[axis].yaxis.set_minor_locator(ticker.AutoMinorLocator())
					# axes_array[axis].yaxis.set_major_locator(ticker.MultipleLocator(5))
					# axes_array[axis].yaxis.set_minor_locator(ticker.MultipleLocator(1))
					# axes_array[axis].yaxis.set_minor_locator(ticker.FixedLocator([0, 10, 0.1]))
					# axes_array[axis].yaxis.set_minor_locator(ticker.LinearLocator(15))
					# axes_array[axis].yaxis.set_major_locator(ticker.IndexLocator(base=.5, offset=.25))
					# axes_array[axis].yaxis.set_yticks([0, 10, 100])
					axes_array[axis].yaxis.set_major_formatter(ScalarFormatter()) # transform logarithm in to continuous ticks
					axis += 1

			# plt.show()
			pdf.savefig(dpi = 75, transparent = True, bbox_inches='tight')
			plt.close()
#################################################################