#!/usr/bin/python3

# Written by Nikhil Goyal, National Institute of Mental Health, 2019-2020
# Purpose: This script does multiple things:
#   1. It creates a .csv file called hcp1200_family_data.csv which is used with CreateKinshipMatrix.R, which generates a .txt file containing a matrix of family relationships
#       1a. This is needed for the 80-20 train-test split testing in the CCA analysis script
#   2. It codifies the fields in the restricted 1200 and unrestricted 1200 files, creating new files called r1200_m.csv and b1200_m.csv, which are to be used in the rest of the analysis pipeline

# usage: python3 sm_file_update.py <../data/b1200.csv> <../data/r1200.csv> <../data/subject_list.txt> <../processed_data/>

import numpy as np
import pandas as pd
from pandas import DataFrame
from numpy import genfromtxt
import os
import sys
from pprint import pprint

cwd = os.getcwd()

b1200_fp = sys.argv[1]
r1200_fp = sys.argv[2]
subjects_fp = sys.argv[3]
out_path = sys.argv[4]


# Import the restricted/behav files for HCP1200
b1200 = pd.read_csv(b1200_fp)
r1200 = pd.read_csv(r1200_fp)

#############
#   PART 2  #
#############
# Now, we need to make a file from which we can use CreateKinshipMatrix.R, needed to create a matrix of twins for proper 80/20 split
# create file which will be used to create a twins list (needed in train-test split)
fields=[
    'Subject',
    'Mother_ID',
    'Father_ID',
    'ZygositySR',
    'ZygosityGT'
]
df=r1200[fields]

# Lastly, drop all but the subjects in iur subject list for HCP1200
subjects = [line.rstrip('.pconn.nii\n') for line in open(subjects_fp)]
family_data = df[df['Subject'].isin(subjects)]

# Save
filename=out_path+'hcp1200_family_data.csv'
out_fp=os.path.join(cwd,filename)
family_data.to_csv(out_fp, index=False)

#############
#   PART 3  #
#############
# Now we will load the r500_m.csv and b500.csv files, and codify their fields as follows:

# Step 1 - Both files:
#   Blank fields --> NaN
#   TRUE/FALSE --> 1/0 (case insensitive),

# Step 2 - Unrestricted file:
#   Gender 
#       F/M --> 1/0
#   Release 
#       Q<1,2,3,...> --> <1,2,3,...>(remove the Q)
#       MEG2  -->  8
#       500  -->  6
#       S900  -->  12
#       S1200 --> 16
#   fMRI_3T_ReconVrs
#       r227 --> 1 
#       anything else --> 0

# Step 3 - Restricted file:
#   Race
#       White                                 -->   0
#       Black or African Am.                  -->   1
#       Asian/Nat. Hawaiian/Othr Pacific Is.  -->   2
#       Am. Indian/Alaskan Nat.               -->   3
#       More than one                         -->   4
#       Unknown or Not Reported               -->   NaN
#   Ethnicity
#       Not Hispanic/Latino         --> 0
#       Hispanic/Latino             --> 1
#       Unknown or Not Reported     --> NaN
#   Color_vision
#       NORMAL  --> 0
#       DEUTAN  --> 1
#       TRITAN  --> 2
#       TRITIAN --> 2
#       PROTAN  --> 3

## STEP 1 - remove blanks, convert true/false ##
# Convert all field values to lowercase (easier to match, case doesn't matter anyway)
b1200 = b1200.astype(str).apply(lambda x: x.str.lower())
r1200 = r1200.astype(str).apply(lambda x: x.str.lower())
# Get rid of any " characters
b1200.apply(lambda x:x.str.replace('"', ""))
r1200.apply(lambda x:x.str.replace('"', ""))
# Remove any leading or trailing whitespaces
b1200 = b1200.applymap(lambda x: x.strip() if isinstance(x, str) else x)
r1200 = r1200.applymap(lambda x: x.strip() if isinstance(x, str) else x)
# Now replace records with ONLY blanks or only spaces with NaN
b1200.replace(r'^\s*$', np.nan, regex=True, inplace=True)
r1200.replace(r'^\s*$', np.nan, regex=True, inplace=True)
# Now convert true/false to 1 or 0
f = lambda x: 1 if x=='true' else (0 if x=='false' else x)
b1200=b1200.applymap(f)
r1200=r1200.applymap(f)


## STEP 2 - Replace values in the UNRESTRICTED file ##
# Gender
b1200['Gender'] = b1200['Gender'].map({'f': 1, 'm': 0}, na_action='ignore')
b1200 = b1200.rename(columns={"Gender":"Sex"})

# Release info
b1200['Release'] = b1200['Release'].str.replace('^q', '',regex=True)  # strip the 'q' from the beginning
# now deal with meg2,s500,s900,s1200
release_dict={'meg2': 8, 's500':6, 's900':12, 's1200':16}
list_of_keys=list(release_dict.keys())
f = lambda x: release_dict[x] if x in list_of_keys else x
b1200['Release'] = b1200['Release'].map(f,na_action='ignore')

# fMRI_3T_ReconVrs
b1200['fMRI_3T_ReconVrs'].replace(to_replace='r227', value=1, regex=False, inplace=True)                 # Replace the EXACT string r227 with 1
b1200['fMRI_3T_ReconVrs'].replace(to_replace=r'^(?!r227$).*', value=0, regex=True, inplace=True)         # Replace anything that doesn't match

# Now, make a copy of this column, but call it 'quarter/release'
b1200['quarter/release'] = b1200['fMRI_3T_ReconVrs']

## STEP 3 - Replace Race, Ethnicity, and Color_Vision fields in restricted file ##
# Race
race={
    'White':0,
    'Black or African Am.':1,
    'Asian/Nat. Hawaiian/Othr Pacific Is.':2,
    'Am. Indian/Alaskan Nat.':3,
    'More than one':4,
    'Unknown or Not Reported':np.nan
    }
# convert dict keys to lowercase
race = {k.lower(): v for k, v in race.items()}
r1200['Race'] = r1200['Race'].map(arg=race, na_action='ignore')

# Ethnicity
ethnicity={
    'Not Hispanic/Latino':0,
    'Hispanic/Latino':1,
    'Unknown or Not Reported':np.nan
}
# convert dict keys to lowercase
ethnicity = {k.lower(): v for k, v in ethnicity.items()}
r1200['Ethnicity'] = r1200['Ethnicity'].map(arg=ethnicity, na_action='ignore')

# Color_vision
vision={
    'NORMAL':0,
    'DEUTAN':1,
    'TRITAN':2,
    'TRITIAN':2,
    'PROTAN':3
}
# convert dict keys to lowercase
vision = {k.lower(): v for k, v in vision.items()}
r1200['Color_Vision'] = r1200['Color_Vision'].map(arg=vision, na_action='ignore')

## STEP 4 - Save these new modified restricted and unrestricted files ##
filename=out_path+'r1200_m.csv'
out_fp=os.path.join(cwd,filename)
r1200.to_csv(out_fp, index=False,na_rep='NaN')

filename=out_path+'b1200_m.csv'
out_fp=os.path.join(cwd,filename)
b1200.to_csv(out_fp, index=False,na_rep='NaN')