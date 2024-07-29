import pandas
import pickle
import re

qc_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/20240501_abcd_fastqc01.txt'
mapping_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/all/ftq_map_mapping.pkl'
output_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/all/20240501_scans.tsv'

# see other good strings to use for exclusions in subset_data.py
exclusions = [
    "Diffusion-FM",
    "Diffusion-FM-AP",
    "Diffusion-FM-PA",
    "fMRI-FM",
    "fMRI-FM-AP",
    "fMRI-FM-PA"
]

# read in the ftq_map files
with open(mapping_file, 'rb') as f:
    mapping = pickle.load(f)

# read in the whole messy abcd_fastqc01.txt file
qc_og = pandas.read_csv(qc_file, low_memory=False, sep='\t')

# drop the messy and unnecessary header line from the dataframe
qc = qc_og.drop(0)

# drop the "ftq_recalled" rows
qc_clean = qc.loc[qc['ftq_recalled'] != '1']
qc_clean = qc_clean.reset_index(drop=True)

# keep only desired columns
intermediary = qc_clean[[
    "subjectkey",
    "src_subject_id",
    "interview_age",
    "sex",
    "visit",
    "file_source",
    "ftq_series_id",
    "abcd_compliant",
    "ftq_complete",
    "ftq_quality",
    "ftq_recalled",
    "ftq_recall_reason",
    "ftq_usable",
    "ftq_notes"
    ]]

mapkeys = list(mapping.keys())
mapvalues = list(mapping.values())

mod_lookup = {
    'bold': 'func',
    'T1w': 'anat',
    'T2w': 'anat',
    'dwi': 'dwi'
}

temp = pandas.DataFrame(data=zip(
            ['/'.join([re.sub(r'_.+', '', x[0]), re.sub(r'sub-.+_ses-(baselineYear1Arm1|[0-9]{1,2}YearFollowUpYArm1)_.+', r'ses-\1', x[0]), mod_lookup[re.sub(r'.+_(.+)\.nii\.gz', r'\1', x[0])], x[0]]) for x in mapvalues],
            [re.sub(r'_.+', '', x[0]) for x in mapvalues],
            [re.sub(r'sub-.+_ses-(baselineYear1Arm1|[0-9]{1,2}YearFollowUpYArm1)_.+', r'ses-\1', x[0]) for x in mapvalues],
            [re.sub(r'.+_(.+)\.nii\.gz', r'\1', x[0]) for x in mapvalues],
            [re.sub(r'.+_(task-[A-z]+)*_.+\.nii\.gz', r'\1', x[0]) if 'task-' in x[0] else 'n/a' for x in mapvalues],
            [re.sub(r'.+_(run-[0-9]+)*_.+\.nii\.gz', r'\1', x[0]) if 'run-' in x[0] else 'n/a' for x in mapvalues],
            mapkeys,
            [re.sub(r'\..+', '', x[0]) for x in mapvalues]
        ), columns=['filename', 'participant_id', 'session_id', 'modality', 'task', 'run', 'ftq_series_id', 'file_prefix'])

merged = pandas.merge(intermediary, temp, on='ftq_series_id', how='outer')

for exclusion in exclusions:
    merged = merged[~merged['ftq_series_id'].str.contains(exclusion)]

final = merged[[
        "filename",
        "participant_id",
        "session_id",
        "modality",
        "task",
        "run",
        "file_prefix",
        "abcd_compliant",
        "ftq_complete",
        "ftq_quality",
        "ftq_recalled",
        "ftq_recall_reason",
        "ftq_usable",
        "ftq_notes",
        "ftq_series_id",
        "file_source"
    ]]

final['session_id'] = pandas.Categorical(final['session_id'], [
    'ses-baselineYear1Arm1',
    'ses-2YearFollowUpYArm1',
    'ses-3YearFollowUpYArm1',
    'ses-4YearFollowUpYArm1',
    'ses-6YearFollowUpYArm1',
    'ses-8YearFollowUpYArm1',
    'ses-10YearFollowUpYArm1'
    ])

final['modality'] = pandas.Categorical(final['modality'], [
    'bold',
    'dwi',
    'T1w',
    'T2w'
    ])

final.sort_values(by=['participant_id', 'session_id', 'modality'], inplace=True)

final.to_csv(output_file, index=False, sep='\t')
