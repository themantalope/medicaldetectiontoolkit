import math
import SimpleITK as sitk
import numpy as np
# from http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/61_Registration_Introduction_Continued.html


def sitk_register_volumes(fixed_image, moving_image, iterations=1000, verbose=True):
    initial_transform = sitk.CenteredTransformInitializer(
        sitk.Cast(fixed_image, moving_image.GetPixelID()),
        moving_image,
        sitk.Euler3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY
    )
    registration = sitk.ImageRegistrationMethod()
    registration.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration.SetMetricSamplingStrategy(registration.RANDOM)
    registration.SetMetricSamplingPercentage(0.01)
    registration.SetInterpolator(sitk.sitkBSpline)
    registration.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=iterations)
    registration.SetOptimizerScalesFromPhysicalShift()
    final_transform = sitk.Euler3DTransform(initial_transform)
    registration.SetInitialTransform(final_transform)

    registration.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    registration.SetSmoothingSigmasPerLevel(smoothingSigmas = [2,1,0])
    registration.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    registration.Execute(sitk.Cast(fixed_image, sitk.sitkFloat32), sitk.Cast(moving_image, sitk.sitkFloat32))
    if verbose:
        print('Optimizer\'s stopping condition, {0}'.format(registration.GetOptimizerStopConditionDescription()))
        print('Final metric value: {0}'.format(registration.GetMetricValue()))
    return moving_image, final_transform



# helper function for image re-sampling
def resample_image(input_image, reference_image, interpolation_method='b-spline', default_intensity=0.0):
    dimension = input_image.GetDimension()
    
    # set up the affine transform between the input and reference
    aff_transform = sitk.AffineTransform(dimension)
    aff_transform.SetMatrix(input_image.GetDirection())
    aff_transform.SetTranslation(np.array(input_image.GetOrigin()) - np.array(reference_image.GetOrigin()))
    
    # modify the transformation to align the centers of the original and reference image instead of their origins.
    centering_transform = sitk.TranslationTransform(dimension)
    
    input_image_center = np.array(
        input_image.TransformContinuousIndexToPhysicalPoint(np.array(input_image.GetSize())/2.0)
    )
    reference_image_center = np.array(
        reference_image.TransformContinuousIndexToPhysicalPoint(np.array(reference_image.GetSize())/2.0)
    )
    
    centering_transform.SetOffset(
        np.array(aff_transform.GetInverse().TransformPoint(input_image_center) - reference_image_center)
    )
    centered_transform = sitk.Transform(aff_transform)
    # centered_transform.AddTransform(centering_transform) 1.x version
    final_transform = sitk.CompositeTransform([centered_transform, centering_transform])
    
    # now we can resample the image
    resampled_image = sitk.Resample(input_image, 
                                    reference_image, 
                                    final_transform, 
                                    interpolation_method, 
                                    default_intensity)
    
    
    return resampled_image


def resample(input_image, new_image_spacing):
    
    img_size = input_image.GetSize()
    img_spacing = input_image.GetSpacing()
    img_origin = input_image.GetOrigin()
    img_dtype = input_image.GetPixelID()
    # img_dimension = input_image.GetDimension()
    img_direction = input_image.GetDirection()

    # get the size of the image
    # new size = (old_spacing/new_spacing) * old_size
    ref_size = [ 
        int(math.floor((osp/nsp) * osi)) for osp, nsp, osi in zip(img_spacing, new_image_spacing, img_size) 
    ]
    # set up the reference image
    ref_img = sitk.Image(ref_size, img_dtype)
    ref_img.SetOrigin(img_origin)
    ref_img.SetSpacing(new_image_spacing)
    ref_img.SetDirection(img_direction)
    
    return resample_image(input_image, ref_img)
    


def sitk_numpy_normalize(sitk_image):
    arr = sitk.GetArrayFromImage(sitk_image)
    mean = arr.mean()
    std = arr.std()
    arr = (arr-mean)/std
    orig_spacing = sitk_image.GetSpacing()
    orig_origin = sitk_image.GetOrigin()
    orig_direction = sitk_image.GetDirection()
    # orig_dtype = sitk_image.GetPixelID()
    # orig_dimension = sitk_image.GetDimension()
    new_img = sitk.GetImageFromArray(arr)
    new_img.SetSpacing(orig_spacing)
    new_img.SetOrigin(orig_origin)
    # new_img.SetPixelID(orig_dtype)
    new_img.SetDirection(orig_direction)

    return new_img