import SimpleITK as sitk
import numpy as np
import os
import resample_itk_image
from itertools import chain
import pandas as pd
import time

# PC
data_folder = r"C:\Users\david\UW\Ryan Kellogg - Kempe files\subjectData"
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

iteration_series = range(0,21,1)
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
    intermed_spacing = 2.0

    image3D_resample = resample_itk_image.resample_img(image3D,[intermed_spacing,intermed_spacing,intermed_spacing])

    print('Resampled Image Spacing ',image3D_resample.GetSpacing())
    print('Resampled Image Size ', image3D_resample.GetSize())

    image3D_resample_seq = resample_itk_image.resample_img(image3D,[min_orig,min_orig,min_orig])

    print('Resampled Image Spacing For Sequential ',image3D_resample_seq.GetSpacing())
    print('Resampled Image Size For Sequential', image3D_resample_seq.GetSize())

    # recast types as necessary for analysis

    recast_image3D = sitk.Cast(sitk.RescaleIntensity(image3D_resample), sitk.sitkUInt8)
    recast32_image3D = sitk.Cast(image3D_resample, sitk.sitkFloat32)
    recast32_image3D_seq = sitk.Cast(image3D_resample_seq, sitk.sitkFloat32)

    #
    initial_seed_point_indexes = np.loadtxt(os.path.join(OUTPUT_DIR,'seeds.csv'),delimiter=',')
    initial_seed_point_indexes = initial_seed_point_indexes.astype(int)
    initial_seed_point_indexes = initial_seed_point_indexes.tolist()


    seg = sitk.Image(image3D_resample.GetSize(), sitk.sitkUInt8)
    seg.CopyInformation(image3D_resample)
    for index in range(1,len(initial_seed_point_indexes)):
        seg[initial_seed_point_indexes[index]] = 1
    seg = sitk.BinaryDilate(seg, 10)

    #stats = sitk.LabelStatisticsImageFilter()
    #stats.Execute(image3D, seg)

    #factor = 3.5
    #lower_threshold = stats.GetMean(1)-factor*stats.GetSigma(1)
    #upper_threshold = stats.GetMean(1)+factor*stats.GetSigma(1)
    #print(stats)

    start_time = time.time()
    # lower and upper intensity thresholds for segmentation
    lower_threshold = 500
    upper_threshold = 1800

    # intialize distance map
    init_ls = sitk.SignedMaurerDistanceMap(seg, insideIsPositive=True, useImageSpacing=True)

    # initialize threshold segmentation level set filter
    lsFilter = sitk.ThresholdSegmentationLevelSetImageFilter()
    lsFilter.SetLowerThreshold(lower_threshold)
    lsFilter.SetUpperThreshold(upper_threshold)
    lsFilter.SetMaximumRMSError(0.007)
    lsFilter.SetNumberOfIterations(1000)
    lsFilter.SetCurvatureScaling(1)
    lsFilter.SetPropagationScaling(1)
    lsFilter.ReverseExpansionDirectionOn()

    #niter = 0
    #for i in range(0,10):
    #    ls = lsFilter.Execute(init_ls,recast32_image3D)
    #    niter += lsFilter.GetNumberOfIterations()
    #    t = "LevelSet after "+str(niter)+" iterations and RMS "+str(lsFilter.GetRMSChange())
    #    fig = myshow(sitk.LabelOverlay(recast_image3D, ls > 0),title=t)

    # run filter
    ls = lsFilter.Execute(init_ls,recast32_image3D)

    # print results
    print(lsFilter)

    elapsed_time = time.time() - start_time
    print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    #myshow(sitk.LabelOverlay(recast_image3D, ls>0))
    # extract >0 values for segmentation
    ls_thresh = ls>0

    image3D_resample_seq_intermed = resample_itk_image.resample_img_set_back_orig(ls_thresh,orig_spacing,orig_size,True)
    image3D_resample_seq = resample_itk_image.resample_img(image3D_resample_seq_intermed,[min_orig,min_orig,min_orig],True)

    start_time = time.time()

    # lower and upper intensity thresholds for segmentation
    lower_threshold = 500
    upper_threshold = 1800

    # intialize distance map
    init_ls_seq = sitk.SignedMaurerDistanceMap(image3D_resample_seq, insideIsPositive=True, useImageSpacing=True)

    # initialize threshold segmentation level set filter
    lsFilterSeq = sitk.ThresholdSegmentationLevelSetImageFilter()
    lsFilterSeq.SetLowerThreshold(lower_threshold)
    lsFilterSeq.SetUpperThreshold(upper_threshold)
    lsFilterSeq.SetMaximumRMSError(0.003)
    lsFilterSeq.SetNumberOfIterations(1000)
    lsFilterSeq.SetCurvatureScaling(1)
    lsFilterSeq.SetPropagationScaling(1)
    lsFilterSeq.ReverseExpansionDirectionOn()

    # run filter
    lsSeq = lsFilterSeq.Execute(init_ls_seq,recast32_image3D_seq)

    # print results
    print(lsFilterSeq)

    elapsed_time = time.time() - start_time
    print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    #myshow(sitk.LabelOverlay(recast_image3D, ls>0))
    # extract >0 values for segmentation
    ls_thresh_seq = lsSeq>0

    ls_thresh_resize = resample_itk_image.resample_img_set_back_orig(ls_thresh_seq,orig_spacing,orig_size,True)
    recast_ls = sitk.Cast(ls_thresh_resize,sitk.sitkUInt8)

    print('Final Image Spacing ',ls_thresh_resize.GetSpacing())
    print('Final Image Size ', ls_thresh_resize.GetSize())

    #shape_stats = sitk.LabelShapeStatisticsImageFilter()
    #shape_stats.ComputeOrientedBoundingBoxOn()
    #shape_stats.Execute(recast_ls)
    #shape_stats.GetPhysicalSize(1)

    # Write the image.
    output_file_name_3D = os.path.join(OUTPUT_DIR, 'segmentation.nii')
    sitk.WriteImage(ls_thresh_resize, output_file_name_3D)
