import pandas as pd

qc_file = '/data/NIMH_scratch/zwallymi/earlea-d2b/fastqc/round_01/abcd_fastqc01_pre_conversion.txt'

# read in the whole messy abcd_fastqc01.txt file
qc_og = pd.read_csv(qc_file, low_memory=False, sep='\t')

# drop the messy and unnecessary header line from the dataframe
qc = qc_og.drop(0)

# drop the "ftq_recalled" rows
qc_clean = qc.loc[qc['ftq_recalled'] != '1']
qc_clean = qc_clean.reset_index(drop=True)

df = qc_clean
ftq_list=[]
for i in df["ftq_series_id"]:
    temp = i.split('_')
    if len(temp) == 4:
        ftq_list.append(temp)
ftq_df = pd.DataFrame(ftq_list, columns=['sub', 'ses', 'series', 'time'])

seslist = []
for i, row in ftq_df.iterrows():
    seslist.append((str(row['sub']), str(row['ses'])))

print(ftq_df["sub"].nunique(), 'unique subjects')
print(len(list(set(seslist))), 'unique sessions')

newdf = pd.DataFrame(list(set(seslist)), columns=['sub', 'ses'])
print(newdf["ses"].value_counts())

sizes = {
    'MID': 227,
    'nBack': 208,
    'SST': 240,
    'rsfMRI': 150,
    'T1': 15,
    'T2': 13,
    'AP': .640,
    'PA': .640,
    'DTI': 114,
    }

#reduces the ftq id to only the scan type
scantype = [x.split("-")[1] for x in ftq_df["series"]]

#determines the file type from each ftq id and adds the corresponding file size from the dictionary
size_in_mb = 0

for i in scantype :
    for key, value in sizes.items() :
        if i.find(key) :
            size_in_mb += value


print(size_in_mb/1024/1024, 'TB of space required')
