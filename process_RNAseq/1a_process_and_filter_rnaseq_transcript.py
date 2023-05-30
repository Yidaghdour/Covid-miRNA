import os

import datatable as dt
import pandas as pd

## CONSTANTS ----------------------------------
tpm_dir = '/scratch/tg1407/COVID_miRNA/data/rnaseq_raw/TPM'

filter_min_tpm = 1
filter_coef = 0.5

# List files
filenames = os.listdir(tpm_dir)
t1_filenames = [s for s in filenames if "T1" in s]

t1_ids = []

for filename in t1_filenames:
    id = filename.split('_')[0]
    t1_ids.append(id)

# For each ID, load and merge
os.chdir(tpm_dir)

for filename in t1_filenames:

    # Identify ID
    id = filename.split('_')[0]

    # Read the filename
    df = pd.read_csv(filename, sep='\t')
    df = df[['transcript_id', 'tpm']]
    df.columns = ['transcript_id', id]

    df = df.groupby('transcript_id').mean().reset_index()

    # Merge
    if filename == t1_filenames[0]:
        rnaseq_df = df
    else:
        rnaseq_df = pd.merge(rnaseq_df, df, on='transcript_id', how='outer')

rnaseq_df

# Filter
for i in range(len(rnaseq_df.index)):
    rnaseq_df.loc[i,'n_above_t'] = sum(rnaseq_df.loc[i,'NGS0709':'NGS0020'] > filter_min_tpm)

rnaseq_filtered_df = rnaseq_df[rnaseq_df['n_above_t'] > filter_coef*len(t1_filenames)]

rnaseq_filtered_df.to_csv('/scratch/tg1407/COVID_miRNA/data/rnaseq_processed/rnaseq_filtered_transcript.csv',
    index=False)
