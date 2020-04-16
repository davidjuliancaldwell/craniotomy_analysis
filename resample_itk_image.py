import SimpleITK as sitk
import numpy as np

def resample_img_set_back_orig(itk_image,orig_spacing,orig_size,is_label=False):

    resample = sitk.ResampleImageFilter()
    resample.SetOutputSpacing(orig_spacing)
    resample.SetSize(orig_size)
    resample.SetOutputDirection(itk_image.GetDirection())
    resample.SetOutputOrigin(itk_image.GetOrigin())
    resample.SetTransform(sitk.Transform())
    resample.SetDefaultPixelValue(itk_image.GetPixelIDValue())


    if is_label:
        resample.SetInterpolator(sitk.sitkNearestNeighbor)
    else:
        resample.SetInterpolator(sitk.sitkBSpline)

    return resample.Execute(itk_image)


# modified from https://gist.github.com/mrajchl/ccbd5ed12eb68e0c1afc5da116af614a
# bspline better ? https://discourse.itk.org/t/which-interpolator-to-resample-image-and-its-segmentation-mask/689/2


def resample_img(itk_image, out_spacing=[2.0, 2.0, 2.0],is_label=False):

    # Resample images to 2mm spacing with SimpleITK
    original_spacing = itk_image.GetSpacing()
    original_size = itk_image.GetSize()

    out_size = [
        int(np.round(original_size[0] * (original_spacing[0] / out_spacing[0]))),
        int(np.round(original_size[1] * (original_spacing[1] / out_spacing[1]))),
        int(np.round(original_size[2] * (original_spacing[2] / out_spacing[2])))]

    resample = sitk.ResampleImageFilter()
    resample.SetOutputSpacing(out_spacing)
    resample.SetSize(out_size)
    resample.SetOutputDirection(itk_image.GetDirection())
    resample.SetOutputOrigin(itk_image.GetOrigin())
    resample.SetTransform(sitk.Transform())
    resample.SetDefaultPixelValue(itk_image.GetPixelIDValue())


    if is_label:
        resample.SetInterpolator(sitk.sitkNearestNeighbor)
    else:
        resample.SetInterpolator(sitk.sitkBSpline)

    return resample.Execute(itk_image)


# modified from https://gist.github.com/mrajchl/ccbd5ed12eb68e0c1afc5da116af614a
# bspline better ? https://discourse.itk.org/t/which-interpolator-to-resample-image-and-its-segmentation-mask/689/2
