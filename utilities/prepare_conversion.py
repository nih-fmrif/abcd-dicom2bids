#! /usr/bin/env python3

# imports
import csv
import datetime
import os
import pandas as pd
import pickle

# from glob import glob
# from collections import Counter
from ftq_map import ftq_map

# set pandas options
pd.set_option('display.max_columns', None)

# inputs
top = '/data/ABCD_DSST/ABCD_BIDS/fast_track/rawdata'
qc_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/20231024_abcd_fastqc01.txt'
output_folder = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/all'
status = os.makedirs(output_folder, exist_ok=True)

# read in the whole messy abcd_fastqc01.txt file
qc_og = pd.read_csv(qc_file, low_memory=False, sep='\t')

# drop the messy and unnecessary header line from the dataframe
qc = qc_og.drop(0)

# drop the "ftq_recalled" rows
qc_clean = qc.loc[qc['ftq_recalled'] != '1']
qc_clean = qc_clean.reset_index(drop=True)

# run the ftq_map function
print(datetime.datetime.now(), 'STARTING FTQ_MAP')
outputs, mapping, errors, have_df, missing = ftq_map(qc_clean['ftq_series_id'], top)
print(datetime.datetime.now(), 'FINISHED FTQ_MAP')

# pickle the outputs
for obj, filename in zip(
    [outputs, mapping, errors, have_df, missing],
    ['ftq_map_outputs.pkl', 'ftq_map_mapping.pkl', 'ftq_map_errors.pkl', 'ftq_map_have_df.pkl', 'ftq_map_missing.pkl']
    ):
    with open(os.path.join(output_folder, filename), 'wb') as f:
        pickle.dump(obj, f)

rm_dirs = {}
not_have = list(have_df[have_df['done'] == 0]['ftq'])
ftq_lst = []
file_lst = []
for i in not_have:
    temp = i.split('_')
    if len(temp) == 4:
        temp.append(i)
        ftq_lst.append(temp)

ftq_df = pd.DataFrame(ftq_lst, columns=['sub', 'ses', 'series', 'time', 'ftq_series_id'])
for subject in ftq_df['sub'].unique():
    sub_df = ftq_df[ftq_df['sub'] == subject]
    temp = []
    for session in sub_df['ses'].unique():
        ses_df = sub_df[sub_df['ses'] == session]
        if os.path.isdir(f'{top}/sub-{subject}/ses-{session}'):
            temp.append(session)
    if temp:
        rm_dirs[subject] = temp

# write out a file of session directories to remove before conversion
with open(os.path.join(output_folder, 'remove_before_conversion.txt'), 'w') as removals:
    for sub in rm_dirs:
        for ses in rm_dirs[sub]:
            removals.write(f'{top}/sub-{sub}/ses-{ses}\n')

errors = []
have = []
not_have = []
drop_idx = []
for index, row in qc_clean.iterrows():
    temp = row['ftq_series_id'].split('_')
    if not len(temp) == 4:
        errors.append(row['ftq_series_id'])
        drop_idx.append(index)
    else:
        if os.path.isdir(f'{top}/sub-{temp[0]}/ses-{temp[1]}'):
            drop_idx.append(index)
            have.append(row['ftq_series_id'])
        else:
            not_have.append(row['ftq_series_id'])

new_qc_clean_drop = qc_clean.drop(drop_idx)
new_qc_clean_drop = new_qc_clean_drop.reset_index(drop=True)
final = pd.concat([qc_og.loc[0:0], new_qc_clean_drop], ignore_index=True)

final.to_csv(os.path.join(output_folder, 'abcd_fastqc01_pre_conversion.txt'), index=False, sep='\t', quoting=csv.QUOTE_ALL)

new_list = list(set(['sub-'+i.replace('_','') for i in list(new_qc_clean_drop['subjectkey'])]))
with open(os.path.join(output_folder, 'full_sublist.txt'), 'w') as file:
    file.write('\n'.join(new_list))
