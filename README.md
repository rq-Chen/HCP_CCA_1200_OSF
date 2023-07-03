# Extending (Smith et al., 2015) to full HCP dataset

This repo is based on the Arxiv preprint [Computationally replicating the Smith et al. (2015) positive-negative mode linking functional connectivity and subject measures](https://www.biorxiv.org/content/10.1101/2020.04.23.058313v2). The authors' github repo is outdated but they provide a link to OSF. This repo is based on the [OSF project of this preprint](https://osf.io/qm49a/).

The files you need (in `data`):

- `RESTRICTED.csv`: the "restricted data" from ConnectomeDB, originally named like `RESTRICTED_yourname_MM_DD_YY_HH_MM_SS.csv`. Requires signing the restricted data use term.
- `UNRESTRICTED.csv`: the "behavioral data" from ConnectomeDB, originally named like `UNRESTRICTED_yourname_MM_DD_YY_HH_MM_SS.csv`.
- `netmats2.txt`: the connectivity matrix from ConnectomDB. It is in the folder `netmats/3T_HCP1200_MSMALL_d200_ts2` after unzipping the file `netmats_3T_HCP1200_MSMAll_ICAd200_ts2.tar.gz` in the `HCP_PTN1200` folder of the downloaded zip file `HCP1200_Parcellation_Timeseries_Netmats.zip`.
- `subject_list.txt` (provided here): ID of the subjects in the netmat file. Originally called `subjectIDs.txt`. Also in the `HCP_PTN1200` folder.
- `column_headers.txt` (provided here): names of the subject measures, obtained from [HCP-CCA website](https://www.fmrib.ox.ac.uk/datasets/HCP-CCA/column_headers.txt) but changing the first variable from 'Subject ID' to 'Subject'.

To run the scripts:

1. Switch to `scripts`.
2. Run the matlab script `summarize_movement.m` to generate a `HCP1200-Motion-Summaries-rfMRI.csv` under the `data` folder. You'll need to modify the 'HCP_dir' in the first line of code.
3. Run the python script `NET.py` by `python3 NET.py ../data/netmats2.txt 200 1003` (netmats location, #ICA components, #subjects) to generate a `NET.txt` under `processed_data`.
4. Run the python script `sm_file_update.py` by `python3 sm_file_update.py ../data/UNRESTRICTED.csv ../data/RESTRICTED.csv ../data/subject_list.txt ../processed_data/` to generate `b1200_m.csv` (modified unrestricted "behavioral" data), `r1200_m.csv` (modified restricted data), `hcp1200_family_data.csv` (family info) under `processed_data`.
5. Run the python script `VARS.py` to generate `VARS.txt` under `processed_data`.
6. Run the matlab script `hcp_1200_cca_final.m` to do the CCA.