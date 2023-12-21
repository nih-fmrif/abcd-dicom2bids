import pandas
import pickle
import re

qc_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/20231024_abcd_fastqc01.txt'
mapping_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/all/ftq_map_mapping.pkl'
output_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/all/20231024_abcd_fastqc01_qc.tsv'

fmaps = [
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

temp = pandas.DataFrame(data=zip(
            [re.sub(r'_.+', '', x[0]) for x in mapvalues],
            [re.sub(r'sub-.+_ses-(baselineYear1Arm1|[246]YearFollowUpYArm1)_.+', r'ses-\1', x[0]) for x in mapvalues],
            [re.sub(r'.+_(.+)\.nii\.gz', r'\1', x[0]) for x in mapvalues],
            [re.sub(r'.+_(task-[A-z]+)*_.+\.nii\.gz', r'\1', x[0]) if 'task-' in x[0] else 'n/a' for x in mapvalues],
            [re.sub(r'.+_(run-[0-9]+)*_.+\.nii\.gz', r'\1', x[0]) if 'run-' in x[0] else 'n/a' for x in mapvalues],
            mapkeys,
            [re.sub(r'\..+', '', x[0]) for x in mapvalues]
        ), columns=['participant_id', 'session_id', 'modality', 'task', 'run', 'ftq_series_id', 'file_prefix'])

merged = pandas.merge(intermediary, temp, on='ftq_series_id', how='outer')

for fmap in fmaps:
    merged = merged[~merged['ftq_series_id'].str.contains(fmap)]

final = merged[[
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

final.to_csv(output_file, index=False, sep='\t')
