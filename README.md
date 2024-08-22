# nirc2_reduction_pipeline

###### Overview ######
This Python objected-oriented project is designed to process astronomical data for a specific target using a sequence of image correction techniques, bad pixel mapping, and distortion corrections. 
The pipeline ultimately produces a final de-rotated and centroid-shifted median image of the target for the purposes of astrometric analysis.

###### Main Features: ######
- Target Setup: Easily define and configure a target system using the Target class.
- Data Cube Creation: Automatically extract frames from FITS files, generating a 3D data cube and computing a 2D median frame.
- Bad Pixel Mapping and Correction: Detects bad pixels using statistical thresholds and applies corrections by replacing their values with the median of surrounding pixels.
- Distortion Correction: Integrates with the RAIN software package to correct image distortions.
- Image De-rotation: Adjusts frames to a common parallactic angle for accurate alignment.
- Centroid Shifting: Aligns images to a common centroid position for further analysis.
- Automated Pipeline Execution: Provides a single command to run the entire sequence from data extraction to generating final processed images.

###### Prerequisites ######
- Python 3.8+
- Required Python Packages:
  - numpy
  - matplotlib
  - astropy
  - photutils
  - scipy
  - skimage


###### Usage ######

Set up your data directory and file ranges for the target system. The directory should contain your FITS files, organized and named sequentially.

Instantiate a Target object with the required parameters and run the entire pipeline:

data_dir_feb = '/path/to/your/data/'
rxs = Target("1rxs1216b", data_dir_feb, [79, 108], [109, 138], [3, 12], [13, 22])
rxs.run_it_all()


###### Outputs ######
The pipeline generates the following output files in a specified directory:

1. Individual frames of the final processed cube as FITS files.
2. A single FITS file containing the final median image.
3. A FITS file containing the final processed data cube.

   
###### Functions Overview ######
Below is a summary of the key functions and their purpose:

get_cube_and_return_median(file_range)
Generates a 3D data cube and computes its 2D median.

create_bad_pixel_map(normalized_flat)
Detects bad pixels using a normalized flat frame.

flag_pixels(index_flags, shape)
Visualizes bad pixels by flagging them on a 2D array.

pixel_filter(array, indices)
Corrects bad pixels by replacing them with the median value of the image.

masked_cube(cube, index_flag)
Applies pixel filtering on a 3D cube of images.

distortion_correction()
Runs the distortion correction code using the RAIN package.

get_final_cube(rain_dir)
Fetches the distortion-corrected images from the RAIN directory.

get_angle_list()
Calculates the angles to de-rotate each frame in a cube.

return_derotated_cube(cube, angle_list)
De-rotates the cube of images to align companions.

shift_centroids(frames, target_location)
Shifts the images based on centroids to a target location.

run_it_all()
Runs all functions sequentially to complete the image processing pipeline.
