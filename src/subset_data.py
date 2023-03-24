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
import subprocess # For calling external programs

from datetime import datetime # For timestamping
from glob import glob # For globbing file names


HERE = os.path.dirname(os.path.realpath(__file__))

__doc__ = """
This command-line tool allows the user to easily subset a directory of ABCC BIDS
input data by datatype (T1, T2, fMRI, field maps).
"""

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-c', '--csv', metavar='CSV' , required=True,
                    help='Input abcd_fastqc01_reformatted.csv file')

parser.add_argument('-i', '--input-dir', metavar='DIRECTORY', required=True,
                    help='Input directory to pull BIDS data from')

parser.add_argument('-o', '--output-dir', metavar='DIRECTORY', required=True,
                    help='Output directory to deposit BIDS hierarchy into')

parser.add_argument('-t', '--types', metavar='TYPE', required=True,
                    nargs='+', choices=DATATYPES.keys(),
                    help='Output directory to deposit BIDS hierarchy into')

args = parser.parse_args()

input_csv = os.path.abspath(args.csv)
input_dir = os.path.abspath(args.input_dir)
output_dir = os.path.abspath(args.output_dir)

# parse datatypes
datatypes = set()
for t in args.types:
    for d in DATATYPES[t]:
        datatypes.add(d)

subsets = list(datatypes)

with open(input_csv, 'r') as f:
    reader = csv.DictReader(f)
    for series in reader:
        # read columns from CSV file
        datatype = series['image_description']
        subject = series['src_subject_id']
        visit = series['EventName']
        age_months = series['SeriesTime']
        sex = series['sex']
        timestamp = series['image_timestamp']
        image_file = series['image_file']
        usable = series['QC']
        compliant = series['ABCD_Compliant']
        complete = series['ftq_complete']
        quality = series['ftq_quality']
        recalled = series['ftq_recalled']
        recall_reason = series['ftq_recall_reason']
        notes = series['ftq_notes']

################################################################################
### LEFT OFF HERE
################################################################################

        # create output TAR.GZ file
        output_basename = '_'.join([nda_bids_subject, bids_session, bids_name, 
                                    'series-' + series_number + '.tar.gz'])
        output = os.path.join(output_dir, 'sourcedata', nda_bids_subject, 
                              bids_session, bids_modality, output_basename)

        # create output directories
        os.makedirs(os.path.dirname(output), exist_ok=True)

        print(datetime.now(), 'Creating TAR.GZ file: ' + output)
        subprocess.run(['tar', '-czf', output, '-C', source, '.'])

# report back to user
print('DICOM to TAR.GZ conversion complete!' + '\n')
print('Output directory: ' + output_dir)

print('Subjects:')
print('\t' + ', '.join(list(set(subjects))) + '\n')

print('Sessions:')
print('\t' + ', '.join(list(set(sessions))) + '\n')

print('Combinations:')
for combo in list(set(combos)):
    print('\t' + ', '.join(combo))
