#!/usr/bin/python3

# Written by Nikhil Goyal, National Institute of Mental Health, 2019-2020
# usage: python3 VARS.py <../data/column_headers.txt> <../data/b1200.csv> <../data/r1200.csv> <../data/subject_list.txt> <../data/HCP1200-Motion-Summaries-rfMRI.csv> <../processed_data/>

# Ex.
# python VARS.py ../data/column_headers.txt ../processed_data/b500_m.csv ../processed_data/r500_m.csv ../data/subject_list.txt ../data/HCP1200-Motion-Summaries-rfMRI.csv ../processed_data/

import numpy as np
import pandas as pd
from pandas import DataFrame
from numpy import genfromtxt
import os
import sys
from pprint import pprint

cwd = os.getcwd()


column_headers_fp = "../data/column_headers.txt"
behavioral_data_fp = "../processed_data/b500_m.csv"
restricted_data_fp = "../processed_data/r500_m.csv"
subject_ids_fp = "../data/subject_list.txt"
motion_data_fp = "../data/HCP1200-Motion-Summaries-rfMRI.csv"
out_path = "../processed_data/"

# column_headers_fp = sys.argv[1]
# behavioral_data_fp = sys.argv[2]
# restricted_data_fp = sys.argv[3]
# subject_ids_fp = sys.argv[4]
# motion_data_fp = sys.argv[5]
# out_path = sys.argv[6]

# get the column headers, and names of subjects
column_headers = [line.rstrip('\n') for line in open(column_headers_fp)]
subjects = [line.rstrip('\n') for line in open(subject_ids_fp)]

# now import "behavioral" and "restricted" datasets into Pandas dataframes
behavioral_data = pd.read_csv(behavioral_data_fp)
restricted_data = pd.read_csv(restricted_data_fp)

# filter the behavioral and restricted datasets to contain only the relevant subject data
behavioral_data = behavioral_data[behavioral_data['Subject'].isin(subjects)]
restricted_data = restricted_data[restricted_data['Subject'].isin(subjects)]

# get the names of column headers
behav_headers=list(behavioral_data.columns.values)
restrict_headers=list(restricted_data.columns.values)

# Make lowercase
column_headers=[element.lower() for element in column_headers]
behav_headers=[element.lower() for element in behav_headers]
restrict_headers=[element.lower() for element in restrict_headers]

# Now let's lets get the column names that are overlapped in each
overlap_in_behav = np.intersect1d(column_headers,behav_headers)
overlap_in_restrict = np.intersect1d(column_headers,restrict_headers)

# Now pull out the columns and their data
# first we will need to convert all the column headers to lowercase
behavioral_data.columns = behavioral_data.columns.str.lower()
restricted_data.columns = restricted_data.columns.str.lower()
behavioral_data_filtered_cols = behavioral_data[overlap_in_behav]
restricted_data_filtered_cols = restricted_data[overlap_in_restrict]

behavioral_data_filtered_cols = behavioral_data[overlap_in_behav]
restricted_data_filtered_cols = restricted_data[overlap_in_restrict]

# Calculate the motion values
motion = pd.read_csv(motion_data_fp)
motion = motion.drop(columns='Scan')
# Each subject has multiple rows (the motion for each of their 4 scans) - average them
motion_sum = motion.groupby('Subject',as_index=False).mean()
# Now drop any subjects not in our subject list
motion_sum = motion_sum[motion_sum['Subject'].isin(subjects)]
motion_sum = motion_sum.rename(columns={"MovementRelativeRMSmean":"rfmri_motion"})

# concat the dataframes
# first reindex all of them to match rfmri_varsqconf
behavioral_data_filtered_cols.index = motion_sum.index
restricted_data_filtered_cols.index = motion_sum.index

frames = [behavioral_data_filtered_cols, motion_sum, restricted_data_filtered_cols]
VARS = pd.concat(frames, axis=1).T.drop_duplicates().T
# VARS = VARS.rename(columns={"subject":"subject id"})
# VARS.reset_index(level=0, inplace=True) 
# get rid of the index for compatibility w/ MATLAB, instead make subject id's in new column called 'subject'
# VARS = VARS.drop(columns='Subject') #drop the subject column so we can reindex w/ column_headers
# VARS.columns.values[0]='Subject' #rename

VARS = VARS.reindex(columns = column_headers)

print("Size of resulting VARS matrix: ", VARS.shape)
filename=out_path+'VARS.txt'
out_fp=os.path.join(cwd,filename)
VARS.to_csv(out_fp, index=False)