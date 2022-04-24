#! /usr/bin/env python3

"""
A tool to subset the ABCD full data set to a smaller subset of datatypes.

Created  2/15/2022 by Eric Earl <eric.earl@nih.gov>
"""

import argparse   # For command line arguments
import csv        # For CSV file handling
import json       # For BIDS sidecar JSON file handling
import os         # For file system operations
import pandas     # For text-file reading and dataframes
import pickle
import subprocess # For calling external programs

from datetime import datetime # For timestamping
from glob import glob # For globbing file names


POSSIBLES = [
    "Diffusion-FM",
    "Diffusion-FM-AP",
    "Diffusion-FM-PA",
    "DTI",
    "fMRI-FM",
    "fMRI-FM-AP",
    "fMRI-FM-PA",
    "MID-fMRI",
    "nBack-fMRI",
    "rsfMRI",
    "SST-fMRI",
    "T1",
    "T1-NORM",
    "T2",
    "T2-NORM"
]

DATATYPES = {
    "anat": [
        "T1",
        "T1-NORM",
        "T2",
        "T2-NORM"
    ],
    "dwi": [
        "Diffusion-FM",
        "Diffusion-FM-AP",
        "Diffusion-FM-PA",
        "DTI",
    ],
    "dwi_fmap": [
        "Diffusion-FM",
        "Diffusion-FM-AP",
        "Diffusion-FM-PA"
    ],
    "func_fmap": [
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA"
    ],
    "fmap": [
        "Diffusion-FM",
        "Diffusion-FM-AP",
        "Diffusion-FM-PA",
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA"
    ],
    "func": [
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA",
        "MID-fMRI",
        "nBack-fMRI",
        "rsfMRI",
        "SST-fMRI"
    ],
    "task-MID": [
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA",
        "MID-fMRI"
    ],
    "task-nback": [
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA",
        "nBack-fMRI"
    ],
    "task-rest": [
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA",
        "rsfMRI"
    ],
    "task-SST": [
        "fMRI-FM",
        "fMRI-FM-AP",
        "fMRI-FM-PA",
        "SST-fMRI"
    ],
    "T1w-asacquired": [
        "T1"
    ],
    "T2w-asacquired": [
        "T2"
    ],
    "T1w-normalized": [
        "T1-NORM"
    ],
    "T2w-normalized": [
        "T2-NORM"
    ]
}


HERE = os.path.dirname(os.path.realpath(__file__))

__doc__ = """
This command-line tool allows the user to easily subset a directory of ABCD-BIDS
input data by datatype (T1, T2, fMRI, field maps, etc).
"""

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-a', '--abcd', metavar='FILE' , required=True,
                    help='Input abcd_fastqc01.txt file')

# parser.add_argument('-i', '--input-dir', metavar='DIRECTORY', required=True,
#                     help='Input directory to pull BIDS data from')

parser.add_argument('-p', '--pickle-file', metavar='FILE', required=True,
                    help='Input path to the "ftq_map_mapping.pkl" file')

parser.add_argument('-o', '--output-dir', metavar='DIRECTORY', required=True,
                    help='Output directory to deposit BIDS hierarchy into')

parser.add_argument('-t', '--types', metavar='TYPE', required=True,
                    nargs='+', choices=DATATYPES.keys(),
                    help="""
                        Space-separated list of data types to subset.  Pick 
                        one or more: """ + ', '.join(DATATYPES.keys())
                    )

parser.add_argument('-g', '--good-qc', action='store_true',
                    help='Only keep data flagged as "good" (ftq_usable==1).')

parser.add_argument('-f', '--intended-for', action='store_true',
                    help='Only keep fmaps with non-empty IntendedFor fields.')

args = parser.parse_args()


# only for testing
# input_dir = os.path.abspath(args.input_dir)
input_dir = os.path.abspath('/data/ABCD_DSST/ABCD_BIDS/fast_track')

rawdata = os.path.join(input_dir, 'rawdata')
sourcedata = os.path.join(input_dir, 'sourcedata')
input_txt = os.path.abspath(args.abcd)
pickle_file = os.path.abspath(args.pickle_file)
output_dir = os.path.abspath(args.output_dir)

with open(pickle_file, "rb") as p:
    ftq_map_mapping = pickle.load(p)

# parse datatypes
datatypes = set()
for t in args.types:
    for d in DATATYPES[t]:
        datatypes.add(d)

subsets = sorted(list(datatypes))

with open(args.abcd, 'r') as f:
    abcd_lines = f.readlines()

for series_idx, series_col in enumerate(abcd_lines[0].split('\t')):
    if 'ftq_series_id' in series_col:
        print('ftq_series_id is in the 1-indexed column: ' + str(series_idx+1))
        break

for recall_idx, recall_col in enumerate(abcd_lines[0].split('\t')):
    if 'ftq_recall_reason' in recall_col:
        print('ftq_recall_reason is in the 1-indexed column: ' + str(recall_idx+1))
        break

for usable_idx, usable_col in enumerate(abcd_lines[0].split('\t')):
    if 'ftq_usable' in usable_col:
        print('ftq_usable is in the 1-indexed column: ' + str(usable_idx+1))
        break

ftq_series_id_dict = {}
for line in abcd_lines[2:]:
    ftq_series_id = line.split('\t')[series_idx].strip('"')
    ftq_recall_reason = line.split('\t')[recall_idx].strip('"')
    ftq_usable = line.split('\t')[usable_idx].strip('"')
    if ftq_recall_reason == '':
        if args.good_qc:
            if ftq_usable == "1":
                pass
            else:
                continue
        ftq_series_id_dict[ftq_series_id] = line

# create output directories
print(datetime.now(), 'Creating subset directory:', output_dir)
os.makedirs(output_dir, exist_ok=True)

# match ftq_series_ids to subsets
print(datetime.now(), 'Matching ftq_series_ids to subsets:', subsets)
ftq_series_ids = []
for subset in subsets:
    print(datetime.now(), '\t', subset)
    for ftq_series_id in ftq_map_mapping.keys():
        if "_ABCD-" + subset + "_" in ftq_series_id:
            ftq_series_ids.append(ftq_series_id)

subset_ftq_series_ids = set(ftq_series_ids)
all_ftq_series_ids = set(ftq_series_id_dict.keys())
final_subset = subset_ftq_series_ids & all_ftq_series_ids

# create subset QC file
subset_qc_file = output_dir + '.abcd_fastqc01.txt'
print(datetime.now(), 'Creating subset QC file:', subset_qc_file)
with open(subset_qc_file, 'w') as f:
    f.write(abcd_lines[0])
    f.write(abcd_lines[1])
    for ftq_series_id in list(final_subset):
        f.write(ftq_series_id_dict[ftq_series_id])

# get the fmaps
print(datetime.now(), 'Collecting relevant field maps, if any')
dwi_fmap_jsons = []
func_fmap_jsons = []
for subset in subsets:
    if 'Diffusion-FM' in subset and dwi_fmap_jsons == []:
        dwi_fmap_jsons = glob(os.path.join(rawdata, 'sub-*', 'ses-*', 'fmap', '*_acq-dwi_*.json'))
    if 'fMRI-FM' in subset and func_fmap_jsons == []:
        func_fmap_jsons = glob(os.path.join(rawdata, 'sub-*', 'ses-*', 'fmap', '*_acq-func_*.json'))

# all fmap JSONs
fmap_jsons = dwi_fmap_jsons + func_fmap_jsons
fmap_jsons_len = len(fmap_jsons)

# in the case of IntendedFor...
final_fmap_jsons = []
if args.intended_for:
    print(datetime.now(), 'Selecting only fmaps with non-empty IntendedFor fields')
    print(datetime.now(), 'This might take a while so please be patient...')
    for i, fmap_json in enumerate(fmap_jsons):
        if i % (round(fmap_jsons_len*0.05)) == 0:
            print(datetime.now(), 'Progress:', str(round(100*i/fmap_jsons_len)), '%')
        with open(fmap_json, 'r') as f:
            fmap_dict = json.load(f)
        if 'IntendedFor' in fmap_dict and not fmap_dict['IntendedFor'] == []:
            final_fmap_jsons.append(fmap_json)

# if no args.intended_for flag was provided then final_fmap_jsons is still empty
if final_fmap_jsons == []:
    # so use all fmap_jsons
    final_fmap_jsons = fmap_jsons

# get the sourcedata
print(datetime.now(), 'Collecting relevant task-based fMRI E-Prime files, if any')
sourcedata_txts = []
for subset in subsets:
    if subset == 'MID-fMRI':
        sourcedata_txts += glob(os.path.join(sourcedata, 'sub-*', 'ses-*', 'func', '*_task-MID_*.txt'))
    if subset == 'nBack-fMRI':
        sourcedata_txts += glob(os.path.join(sourcedata, 'sub-*', 'ses-*', 'func', '*_task-nback_*.txt'))
    if subset == 'SST-fMRI':
        sourcedata_txts += glob(os.path.join(sourcedata, 'sub-*', 'ses-*', 'func', '*_task-SST_*.txt'))

# symlink the ftq_series_id mapped files
print(datetime.now(), 'Symbolically linking func, dwi, and anat files, if any')
for ftq_series_id in list(final_subset):
    if '-fMRI_' in ftq_series_id or '_ABCD-rsfMRI_' in ftq_series_id:
        modality = 'func'
    elif '_ABCD-DTI_' in ftq_series_id:
        modality = 'dwi'
    elif '_ABCD-T1' in ftq_series_id or '_ABCD-T2' in ftq_series_id:
        modality = 'anat'

    mapped_file_list = ftq_map_mapping[ftq_series_id]
    subses_underscore_split = mapped_file_list[0].split('_')
    subject = subses_underscore_split[0]
    session = subses_underscore_split[1]
    for mapped_file in mapped_file_list:
        mapped_path = os.path.join(rawdata, subject, session, modality, mapped_file)
        mapped_file_relpath = os.path.relpath(mapped_path, input_dir)
        output_path = os.path.join(output_dir, mapped_file_relpath)
        output_subdir = os.path.dirname(output_path)
        os.makedirs(output_subdir, exist_ok=True)
        os.symlink(mapped_path, output_path)

# symlink the fmaps, if any
print(datetime.now(), 'Symbolically linking fmap files, if any')
for fmap_json in final_fmap_jsons:
    fmap_nifti = fmap_json.replace('.json', '.nii.gz')
    fmap_json_relpath = os.path.relpath(fmap_json, input_dir)
    output_json = os.path.join(output_dir, fmap_json_relpath)
    output_nifti = os.path.join(output_dir, fmap_json_relpath.replace('.json', '.nii.gz'))
    output_subdir = os.path.dirname(output_json)
    os.makedirs(output_subdir, exist_ok=True)
    os.symlink(fmap_json, output_json)
    os.symlink(fmap_nifti, output_nifti)
    if '_acq-dwi_' in fmap_json:
        fmap_bval = fmap_json.replace('.json', '.bval')
        fmap_bvec = fmap_json.replace('.json', '.bvec')
        output_bval = os.path.join(output_dir, fmap_json_relpath.replace('.json', '.bval'))
        output_bvec = os.path.join(output_dir, fmap_json_relpath.replace('.json', '.bvec'))
        os.symlink(fmap_bval, output_bval)
        os.symlink(fmap_bvec, output_bvec)

# symlink the sourcedata, if any
print(datetime.now(), 'Symbolically linking task-based fMRI E-Prime files, if any')
for sourcedata_txt in sourcedata_txts:
    sourcedata_txt_relpath = os.path.relpath(sourcedata_txt, input_dir)
    output_txt = os.path.join(output_dir, sourcedata_txt_relpath)
    output_subdir = os.path.dirname(output_txt)
    os.makedirs(output_subdir, exist_ok=True)
    os.symlink(sourcedata_txt, output_txt)

print(datetime.now(), 'All done!')
