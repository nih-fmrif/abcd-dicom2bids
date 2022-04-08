#! /usr/bin/env python3

"""
A tool to subset the ABCC data to a smaller set of datatypes.

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
This command-line tool allows the user to easily subset a directory of ABCC BIDS
input data by datatype (T1, T2, fMRI, field maps).
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
                    help='Space-separated list of data types to subset.  Pick one or more: ' + ', '.join(DATATYPES.keys()))

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

for index, column in enumerate(abcd_lines[0].split('\t')):
    if 'ftq_series_id' in column:
        print('ftq_series_id is in the 1-indexed column: ' + str(index+1))
        break

ftq_series_id_dict = {}
for line in abcd_lines[2:]:
    ftq_series_id = line.split('\t')[index].strip('"')
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

# abcd_df = pandas.read_csv(args.abcd, sep='\t', header=0, skiprows=[1])

# with open(input_csv, 'r') as f:
#     reader = csv.DictReader(f)
#     for series in reader:
#         # read columns from CSV file
#         datatype = series['image_description']
#         subject = series['src_subject_id']
#         visit = series['EventName']
#         age_months = series['SeriesTime']
#         sex = series['sex']
#         timestamp = series['image_timestamp']
#         image_file = series['image_file']
#         usable = series['QC']
#         compliant = series['ABCD_Compliant']
#         complete = series['ftq_complete']
#         quality = series['ftq_quality']
#         recalled = series['ftq_recalled']
#         recall_reason = series['ftq_recall_reason']
#         notes = series['ftq_notes']
