import SimpleITK as sitk
import numpy as np
import os
import resample_itk_image
from itertools import chain
import pandas as pd
import time

# PC
data_folder = r"C:\Users\david\UW\Ryan Kellogg - Kempe files\subject_data"
data_folder = os.path.abspath(os.path.expanduser(os.path.expandvars(data_folder)))

# create output directories if they don't exist

#for root, dirs, files in os.walk(data_folder, topdown=True):
    # walking tips from https://stackoverflow.com/questions/42720627/python-os-walk-to-certain-level
#    if root[len(data_folder):].count(os.sep) == 2:
#        if not (os.path.exists(os.path.join(root,'output'))):
#            os.mkdir(os.path.join(root,'output'))

# generate structure of files

data_subj_list = list()
output_dir_list = list()
for data_subj, dirs, files in os.walk(data_folder, topdown=True):
    # walking tips from https://stackoverflow.com/questions/42720627/python-os-walk-to-certain-level
    if (data_subj[len(data_folder):].count(os.sep) == 3) and not 'output' in data_subj[len(data_folder):]:
       # print(data_subj)
        super_folder = data_subj.rsplit('\\',1)
        OUTPUT_DIR = os.path.join(super_folder[0],'output')

        data_subj_list.append(data_subj)
        output_dir_list.append(OUTPUT_DIR)

length_data = len(data_subj_list)

iteration_series = range(0,length_data,1)
#iteration_series = [0]

for select in iteration_series:
    # Read the original series. First obtain the series file names using the
    # image series reader.
    data_subj = data_subj_list[select]
    OUTPUT_DIR = output_dir_list[select]

    print('Subject folder = ',data_subj)
    series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(data_subj)
    if not series_IDs:
        print("ERROR: given directory \""+data_folder+"\" does not contain a DICOM series.")
    print('File # ',select)
    print('\n')
