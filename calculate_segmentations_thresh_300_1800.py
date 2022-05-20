import SimpleITK as sitk
import numpy as np
import os
import resample_itk_image
from itertools import chain
import pandas as pd
import time

lower_threshold = 300
upper_threshold = 1800

# PC
data_folder = r"C:\Users\david\UW\Ryan Kellogg - Kempe files\subjectData"
data_folder = os.path.abspath(os.path.expanduser(os.path.expandvars(data_folder)))

# create output directories if they don't exist

#for root, dirs, files in os.walk(data_folder, topdown=True):
   # walking tips from https://stackoverflow.com/questions/42720627/python-os-walk-to-certain-level
#   if root[len(data_folder):].count(os.sep) == 2:
#       if not (os.path.exists(os.path.join(root,'output'))):
#           os.mkdir(os.path.join(root,'output'))

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

iteration_series = range(3,5,1)
#iteration_series  =[3]

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

    series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(data_subj, series_IDs[0])

    series_reader = sitk.ImageSeriesReader()
    series_reader.SetFileNames(series_file_names)

    # Configure the reader to load all of the DICOM tags (public+private):
    # By default tags are not loaded (saves time).
    # By default if tags are loaded, the private tags are not loaded.
    # We explicitly configure the reader to load tags, including the
    # private ones.
    series_reader.MetaDataDictionaryArrayUpdateOn()
    #series_reader.LoadPrivateTagsOn()
    image3D = series_reader.Execute()

    reader = sitk.ImageFileReader()
    reader.SetFileName(series_file_names[0])
    #reader.LoadPrivateTagsOn();

    reader.ReadImageInformation();

    print('Original Image Spacing ',image3D.GetSpacing())
    print('Original Image Size ', image3D.GetSize())

    orig_size = image3D.GetSize()
    orig_spacing = image3D.GetSpacing()
    min_orig = min(orig_spacing)
    max_orig = max(orig_spacing)

    if max_orig > 50:
        print('Error: File ',select,' may have extra dicom file, move -0001 to another folder, try again','\n')
        exit()

    image3D_resample_seq = resample_itk_image.resample_img(image3D,[min_orig,min_orig,min_orig])

    print('Resampled Image Spacing For Sequential ',image3D_resample_seq.GetSpacing())
    print('Resampled Image Size For Sequential', image3D_resample_seq.GetSize())

    # recast types as necessary for analysis

    start_time = time.time()

    seg = sitk.BinaryThreshold(image3D_resample_seq, lowerThreshold=lower_threshold, upperThreshold=upper_threshold, insideValue=1, outsideValue=0)

    # binary closing to smooth out skull
    if min_orig <0.6:
        vectorRadius=(3,3,3)
    else:
        vectorRadius=(2,2,2)

    print('Morphological Closing radius = ',vectorRadius)

    kernel=sitk.sitkBall
    seg_resized_clean = sitk.BinaryMorphologicalClosing(seg, vectorRadius, kernel)

    # find biggest connected object, which is skull, get rid of rest.
    connected = sitk.ConnectedComponent(seg_resized_clean)
    connected = sitk.RelabelComponent(connected,sortByObjectSize=True)
    connected_select = connected == 1

    final_cleaned = resample_itk_image.resample_img_set_back_orig(connected_select,orig_spacing,orig_size,True)

    elapsed_time = time.time() - start_time
    print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    print('\n')

    # Write the image.
    output_file_name_3D = os.path.join(OUTPUT_DIR, 'segmentation_thresh_300_1800.nii')
    sitk.WriteImage(final_cleaned, output_file_name_3D)
