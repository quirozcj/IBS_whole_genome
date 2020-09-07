import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from numpy import *

 
sample_file = (sys.argv[1])
positive_control_haplotypes_file = (sys.argv[2])
wind_size = (sys.argv[3])
components_no = int((sys.argv[4]))
sample_ref = (sys.argv[5])
sample_query = (sys.argv[6])
assembly_suffix = (sys.argv[7])
FILTER=int((sys.argv[8]))

sys.stdout = open(sys.argv[9],'w')
sys.stdout2 = open(sys.argv[10],'w')
sys.stdout3 = open(sys.argv[11],'w')
# sys.stdout4 = open(sys.argv[12],'w')

## output file design for precision-recall summary
## NOTE: PR = short for precision-recal
## NOTE: df = short for DataFrame
PR_summary = []
header ='filter''\t' \
            'sum_TP''\t' \
            'sum_FP''\t' \
            'sum_TN''\t' \
            'sum_FN''\t' \
            'precision''\t' \
            'recall''\t' \
            'F1''\t' \
            'accuracy''\t' \
            'specificity''\n'
PR_summary.append(header)
## generate file for precision-recall summary data

#when combined files
file_output = 'window size: ' + str(wind_size) + '\n'

# indivudual files
# PR_out_file_PR_summary_ind = open('precision_recall_summary_whole_genome_{0}_{1}_{2}_{3}_kb_wind_gauss_comp{4}.txt'.format((sample_ref), (sample_query), (assembly_suffix), (wind_size), (components_no)), 'w')
# PR_out_file_PR_summary_ind

# PR_out_file_PR_summary = open('precision_recall_summary_whole_genome_{0}_{1}_{2}_kb_wind_gauss_comp{3}.txt'.format((sample_ref), (sample_query), (assembly_suffix), (components_no)), 'a')
# PR_out_file_PR_summary.write(file_output) # header for combined files

## filter for high haplotig-count
# filter_range = range(200, 600, 200)                                                                            				# filter used end not inlcusive
# for FILTER in filter_range:
    ## create files for gaussian model and precision-recall output (one per range)
    # NOTE: I disabled create the original gauss_mx file since we are using only the normalized file
    # gauss_PR_out_file = open('whole_genome_haplo-tig_count_{0}_{1}_{2}_{3}_kb_wind_filter_{4}_gauss_comp{5}.txt'.format((sample_ref), (sample_query), (assembly_suffix), (wind_size), (FILTER), (components_no)), 'w')
    # PR_out_file = open('precision_recall_whole_genome_haplo-tig_count_{0}_{1}_{2}_{3}_kb_wind_filter_{4}_gauss_comp{5}.txt'.format((sample_ref), (sample_query), (assembly_suffix), (wind_size), (FILTER), (components_no)), 'w')


    ############## generate gaussian_mx grouping model and file ###################

    ## get the data ready
input_file = pd.read_csv(sample_file, delimiter='\t')
filtered_df_gauss = input_file[input_file['haplotig_count']<= FILTER]                                                  # filter per count
filtered_df_gauss = pd.DataFrame(filtered_df_gauss[['chr','start','end','haplotig_count']])                             # get needed columns
filtered_df_gauss.reset_index(drop=True, inplace=True)                                                            # reset index to match calculations
haplotig_count_array = np.array(filtered_df_gauss['haplotig_count'])                                              # dataframe filtered into array
log_haplotig_count_array = np.log(haplotig_count_array, where=(haplotig_count_array!=0))
log_haplotig_count_array = np.nan_to_num(log_haplotig_count_array, nan=0.0)
    # where_are_NaNs = isnan(log_haplotig_count_array)
    # log_haplotig_count_array[where_are_NaNs] = 0                                                       # transform the data into log scale
log_haplotig_count_array = log_haplotig_count_array.reshape((len(log_haplotig_count_array),1))       # make array in two dimensions
log_haplotig_count_array = np.nan_to_num(log_haplotig_count_array, nan=0.0)
    # where_are_NaNs = isnan(log_haplotig_count_array)
    # log_haplotig_count_array[where_are_NaNs] = 0                         

    ### Model
model = GaussianMixture(n_components=components_no, covariance_type='full')                                   # build model, covariance default is 'full'
model.fit(log_haplotig_count_array)                                                               # Estimate model parameters with the EM algorithm (trainning the model)
log_haplotig_count_array_fit = model.fit(log_haplotig_count_array)
array_predicted_model_log = model.predict(log_haplotig_count_array)                               # Predict the labels for the data samples in X using trained model (I am using the same trainnin to predict)

array_predicted_model_log = array_predicted_model_log.reshape((len(array_predicted_model_log),1)) # two dimensions array
array_predicted_model_log = pd.DataFrame({'gauss_mx': array_predicted_model_log[:, 0]})           # dataframe form array
filtered_df_gauss['gauss_mx'] = array_predicted_model_log['gauss_mx']                                   # combine dataframe 
    
    ### write gaussian model file. NOW is disabled
    # filtered_df_gauss.to_csv(gauss_PR_out_file, index=False,  sep='\t')                # write output to a file
filtered_df_gauss.to_csv(sys.stdout, index=False,  sep='\t')              # write output to the terminal




    ############## Positive Control #################
    ## open positive control haplotypes file (Brinton et al 2020)
    ## get the data ready
input_file_1 = pd.read_csv(positive_control_haplotypes_file, delimiter='\t')
filtered_df = input_file_1[(input_file_1['ref'] == sample_ref) & (input_file_1['query'] == sample_query)]                         # filter per count
filtered_df = pd.DataFrame(filtered_df[['ref','query','chrom','ref_start','ref_end']]).sort_values(['chrom','ref_start'])   # get needed columns
filtered_df.reset_index(drop=True, inplace=True)                                                           				# reset index to match calculations

    ## Use kmer-haplotigs gaussian prediciton DataFrame(filtered_df_gauss) to get precision-recall
    ## normalize IBS from gaussian groups into 1s (IBS) and 0s (NON-IBS), (dummy the data)
get_IBS = filtered_df_gauss.groupby(by=['gauss_mx'])
key_group=[]
value_group=[]
for k, group in get_IBS:
    key_group.append(k)
    value_group.append(group['haplotig_count'].median())
dic_grupby = dict(zip(key_group, value_group))
IBS_group = min(dic_grupby, key=dic_grupby.get)
filtered_df_gauss['gauss_mx'] = np.where(filtered_df_gauss['gauss_mx'] == IBS_group, 1, 0)  # condition where "gauss_mx" values = to IBS_group, "1" (IBS) will be assigned,
                                                                                                # for all the reamaining (NON-IBS) "0" will be given

    ## merge data frames information from gauss_mx and align-IBS (Brinton et all 2020) DataFrames)
    ## get IBS only
match_chrom = pd.merge(filtered_df_gauss, filtered_df, left_on='chr', right_on='chrom', how='right')
match_chrom = match_chrom.loc[(match_chrom['start'] >= match_chrom['ref_start']) & (match_chrom['end'] < match_chrom['ref_end'])]
match_chrom['align_IBS'] = pd.Series(1, index=match_chrom.index)
match_chrom.reset_index(drop=True, inplace=True)
    # use the IBS region obtained to get all IBS and NON-IBS across the genome
merge_non_match = pd.merge(filtered_df_gauss, match_chrom, left_on=['chr', 'start'], right_on=['chr', 'start'], how='left', suffixes=('_left', '_right'))

    ## get categorical values for precision-recall analysis
merge_non_match['TP'] = ((merge_non_match['gauss_mx_left']==1) & (merge_non_match['align_IBS']==1)).astype(int)
merge_non_match['FP'] = ((merge_non_match['gauss_mx_left']==1) & (merge_non_match['align_IBS']!=1)).astype(int)
merge_non_match['TN'] = ((merge_non_match['gauss_mx_left']!=1) & (merge_non_match['align_IBS']!=1)).astype(int)
merge_non_match['FN'] = ((merge_non_match['gauss_mx_left']!=1) & (merge_non_match['align_IBS']==1)).astype(int)
merge_non_match['align_IBS'] = ((merge_non_match['align_IBS']==1)).astype(int) 				# transform IBS from align into "1" and "0"
merge_non_match['ref'] = (filtered_df['ref'][1])        									# get ref and query sample name
merge_non_match['query'] = (filtered_df['query'][1])    
    # pd.set_option('display.max_columns', None) # print full columns

    ## extract useful data and rename columns
PR_df = pd.DataFrame(merge_non_match[['ref', 'query', 'chr', 'start', 'end_left', \
                                                        'haplotig_count_left','gauss_mx_left', \
                                                        'align_IBS', \
                                                        'TP',  'FP',  'TN','FN' ]])

PR_df.rename(columns={'haplotig_count_left':'haplotig_count',
                                        'gauss_mx_left':'gauss_mx_IBS', \
                                        'end_left':'end'}, inplace=True)

    ## write into ouput file with categorical data as "1" (IBS) and "0" as (NON-IBS)
    # PR_df.to_csv(PR_out_file, index=False,  sep='\t')                	# write output to a file
    # PR_df.to_csv(sys.stdout, index=False,  sep='\t')              # write output to the terminal
PR_df.to_csv(sys.stdout2, index=False,  sep='\t')




    ############# precision-recall summary #############

## filter and remove "chrUn" for downstream analysis
PR_df = PR_df.loc[PR_df['chr'] != 'chrUn']

    ## get recall statistics using the "PR_df" DataFrame
sum_TP = PR_df['TP'].sum()
sum_FP = PR_df['FP'].sum()
sum_TN = PR_df['TN'].sum()
sum_FN = PR_df['FN'].sum()

    # append results to a list, to create DataFrame later on
PR_summary.append(str(FILTER)+'\t'+ str(sum_TP) +'\t'+ str(sum_FP) +'\t'+ str(sum_TN) +'\t'+ str(sum_FN)+'\n')


# create dataframe from a list to calculate precision-recall
PR_list=[]
for x in PR_summary:
	PR_list.append(x.rstrip().split('\t'))
df_PR_summary = pd.DataFrame(PR_list)
new_header = df_PR_summary.iloc[0]
df_PR_summary = (df_PR_summary[1:])
df_PR_summary.columns = new_header

# transform columns df to numeric (int)
df_PR_summary['sum_TP'] = pd.to_numeric(df_PR_summary['sum_TP'])
df_PR_summary['sum_FP'] = pd.to_numeric(df_PR_summary['sum_FP'])
df_PR_summary['sum_TN'] = pd.to_numeric(df_PR_summary['sum_TN'])
df_PR_summary['sum_FN'] = pd.to_numeric(df_PR_summary['sum_FN'])

# Precision: TP/TP+FP
df_PR_summary['precision'] = \
(df_PR_summary['sum_TP']/ \
(df_PR_summary['sum_TP']+ \
df_PR_summary['sum_FP'])).round(3)

# Recall: TP/TP+FN
df_PR_summary['recall'] = \
(df_PR_summary['sum_TP']/ \
(df_PR_summary['sum_TP']+ \
df_PR_summary['sum_FN'])).round(3)

# F1: 2*((precision * recall)/(precision + recall))
df_PR_summary['F1'] = \
(2*((df_PR_summary['precision']* \
df_PR_summary['recall'])/ \
(df_PR_summary['precision']+ \
df_PR_summary['recall']))).round(3)

# Accuracy: TP+TN/TP+FN+TN+FP
df_PR_summary['accuracy'] = \
((df_PR_summary['sum_TP']+ \
df_PR_summary['sum_TN'])/ \
(df_PR_summary['sum_TP']+ \
df_PR_summary['sum_FN']+ \
df_PR_summary['sum_TN']+ \
df_PR_summary['sum_FP'])).round(3)

# Specificity: TN/TN+FP
df_PR_summary['specificity'] = \
(df_PR_summary['sum_TN']/ \
(df_PR_summary['sum_TN']+ \
df_PR_summary['sum_FP'])).round(3)

# write summary to a file
df_PR_summary.to_csv(sys.stdout3, index=False,  sep='\t')
# df_PR_summary.to_csv(PR_out_file_PR_summary, index=False,  sep='\t')
# PR_out_file_PR_summary.write('\n')
# df_PR_summary.to_csv(sys.stdout, index=False,  sep='\t')