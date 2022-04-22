#! /usr/bin/env python3

"""
A tool to subset the ABCD full data set to a smaller subset of datatypes.

Created  2/15/2022 by Eric Earl <eric.earl@nih.gov>
"""

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

import argparse   # For command line arguments
import csv        # For CSV file handling
import os         # For file system operations
import pandas     # For text-file reading and dataframes
import pickle
import subprocess # For calling external programs

from datetime import datetime # For timestamping
from glob import glob # For globbing file names


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

args = parser.parse_args()

input_txt = os.path.abspath(args.abcd)
# input_dir = os.path.abspath(args.input_dir)
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
