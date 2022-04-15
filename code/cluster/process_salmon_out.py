import pandas as pd
import numpy as np
from dppd import dppd
import glob
import os
import sys

'''
description
'''

print('- initialize tpm and raw count matricies')
dataset = sys.argv[1]
print('-- dataset name: ' + dataset)
function = sys.argv[2]
print('-- function name: ' + function)
a_quant_file = glob.glob('processed/salmon_out/' + dataset + '/' + '/*/quant.sf')[0]
quant_file = pd.read_csv(a_quant_file, sep = '\t', header = 0)
tpm = pd.DataFrame(quant_file[['Name', 'Length', 'EffectiveLength']])
counts = pd.DataFrame(quant_file[['Name', 'Length', 'EffectiveLength']])

print('- iterating over quant files and fill tpm and count matricies')
# filenames = [i for i in glob.glob('processed/salmon_out/*/quant.sf')]
# quants = [pd.read_csv(file, sep = "\t") for file in filenames]

for file in glob.glob('processed/salmon_out/' + dataset + '/*/quant.sf'):
	
	# get current run and quant info
	current_run = file.split('/')[3]
	current_quant = pd.read_csv(file, sep = "\t")
	
	# add current info to master matricies
	tpm[current_run] = current_quant['TPM']
	counts[current_run] = current_quant['NumReads']

print('- saving TPM and count matricies')

print('-- TPM file path:  ' + 'results/salmon_processed/' + dataset + '/' + function + '_' + 'tpm.tsv')
tpm.to_csv('results/salmon_processed/' + dataset + '/' + function + '_' + 'tpm.tsv', sep = '\t', index = False)
print('-- counts file path:  ' + 'results/salmon_processed/' + dataset + '/' + function + '_' + 'counts.tsv', sep = '\t')
counts.to_csv('results/salmon_processed/' + dataset + '/' + function + '_' + 'counts.tsv', sep = '\t', index = False)
