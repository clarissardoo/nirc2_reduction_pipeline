import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import sys
import photutils
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)
import scipy.ndimage

from skimage.transform import rotate as skirotate


class Target:
    '''
    SET UP A TARGET CLASS, WHERE TARGET IS SPECIFIC SYSTEM. THE TARGET REQUIRES:

    name: string with target name
    data_dir: string of directory where data is located
    file_range: array with the format [initial, final] for file numbers with target frames
    sky_filerange: same as file_range, but for sky frames
    flaton_filerange: same as file_range, but for flats on frames
    flatoff_filerange: ssame as file_range, but for flats off frames

    an example of a target setup is: rxs = Target("1rxs1216b", data_dir_feb, [79, 108], [109, 138], [3, 12], [13, 22])

    '''
    def __init__(self, name, data_dir, file_range, sky_filerange, flaton_filerange, flatoff_filerange):
        self.name = name
        self.data_dir = data_dir
        self.file_range = file_range
        self.sky_filerange = sky_filerange
        self.flaton_filerange = flaton_filerange
        self.flatoff_filerange = flatoff_filerange


    def get_cube_and_return_median(self, file_range):
        '''
        This function obtains a list of frames and outputs a cube and its 2d-median. 
        The function already divides the frames by number of coadds

        Args: 
        file_range: array [initial, final] where initial and final are ints of file numbers with target frames

        Returns:

        cube - 3d numpy array with shape [frame number, pixel x, pixel y]
        cube_median - median 2d numpy array with shape [pixel x, pixel y]

        '''
        data_files = np.array(sorted(os.listdir(self.data_dir)))
        cube = []
        for i in data_files[file_range[0]-1:file_range[-1]]:
            #print(file_range)
            frame = fits.getdata(self.data_dir + i)
            h = fits.open(self.data_dir +i)
            coadds = h[0].header['COADDS']
            cube.append(frame/coadds)
        cube = np.array(cube)
        cube_median = np.median(cube, axis = 0)

        return cube, cube_median

    def create_bad_pixel_map(self,normalized_flat):
        '''
        This function obtains bad pixel locations provided a normalized flat frame. 

        Args: 

        normalized_flat: 2d numpy array of normalized flat frame

        Returns:

        index_flags: numpy array with tuples corresponding to locations of bad pixels

        '''

        stdev = np.std(normalized_flat)
        mean = np.mean(normalized_flat)

        five_sig_lower = mean - 5*stdev

        index_flags = []
        for i, j in np.argwhere(normalized_flat < five_sig_lower): 
            index_flags.append((i, j))

        five_sig_higher = mean + 5*stdev

        for i, j in np.argwhere(normalized_flat > five_sig_higher): 
            index_flags.append((i, j))
        return index_flags

    def flag_pixels(self,index_flags,shape):

        '''
        This function flags bad pixels on a 2d numpy array. 
        Mainly a visualization function; not used in the pipeline

        Args: 

        index_flags: output from create_bad_pixel_map (numpy array with tuples corresponding to location of bad pixels)
        shape: shape of desired image. Formatting is a tuple (x shape, y shape)

        Returns:

        array: 2d numpy array where good pixels are 0's and bad pixels are 1's

        '''
        array = np.zeros(shape, dtype=int)
        for index in index_flags:
            row, col = index
            array[row, col] = 1
        return array

    def pixel_filter(self,array, indices):

        '''
        This filters bad pixel values by replacing their count with median value of image

        Args: 

        array: 2d-array mapping where pixel filter needs to be applied
        indices: output from create_bad_pixel_map (numpy array with tuples corresponding to location of bad pixels) 

        Returns:

        fixed_array: 2d pixel-corrected array 

        '''
        for index in indices:
            row, col = index
            median = np.median(array)
            fixed_array[row, col] = median
        return fixed_array

    def masked_cube(self,cube, index_flag):

        '''
        This filters bad pixel values by replacing their count with median value of image on a CUBE

        Args: 

        cube: 3d-array of cube where correction needs to be applied
        indices: output from create_bad_pixel_map (numpy array with tuples corresponding to location of bad pixels) 

        Returns:

        fixed_array: 3d pixel-corrected cube

        '''
        cube_masked = []
        for i in range(cube.shape[0]):
            target_corrected = self.pixel_filter(cube[i],index_flag)
            cube_masked.append(target_corrected)
        return cube_masked

    def distortion_correction(self):
        '''
        This runs distortion correction code using rain

        Args: 

        No args

        Returns:

        rain_dir: rain directory path where undistorted pictures are located

        '''
        path = '/home/cdoo/rain_new/'
        os.chdir(path)
        pythoncode = os.system('mkdir ' + self.name)
        pythoncode2 = os.system('cp -r ' + '/home/cdoo/nirc2_data/' +str(self.name) + ' ' + str(path) + str(self.name))


        new_path = path + '/' + str(self.name)

        code_path = '/home/cdoo/rain_new/rain/'
        os.chdir(code_path)

        rain_code = os.system('python3 rain_keck_data.py ' + path + str(self.name) + '/' + str(self.name) + ' ' + str(self.name) )
        print('python3 rain_keck_data.py ' + path + '/' + str(self.name) + ' ' + str(self.name)  )

        rain_dir = path+self.name+ '/' + self.name


        return rain_dir


    def get_final_cube(self, rain_dir):
        '''
        This gets distortion corrected images from rain directory

        Args: 

        rain_dir: directory path output from distortion_correction (string format)

        Returns:

        rain_dir: rain directory path where undistorted pictures are located

        '''
        rain_dir = rain_dir
        rain_files = np.array(sorted(os.listdir(rain_dir)))
        if rain_files[0] == '.DS_Store':
            rain_files = rain_files [1:]

        #run code, then get data here
        final_rain_files = []
        for rain_file in rain_files:
            if "rain_" in rain_file:
                final_rain_files.append(rain_file)

        final_cube = []
        for i in final_rain_files:
            #print(i)
            frame = fits.getdata(rain_dir + '/' + i, ignore_missing_simple=True)
            print(frame)
            final_cube.append(frame)

        final_cube = np.array(final_cube)

        return final_cube



    def get_angle_list(self):

        '''
        This gets angle that pictures in a cube need to be de-rotated by

        Args: 

        No args

        Returns:

        angle_list: python list with angles to distort each image in a cube

        '''
        data_files = np.array(sorted(os.listdir(self.data_dir)))
        angle_list = []
        zp_offset = -0.262
        for i in data_files[self.file_range[0]-1:self.file_range[1]]:
            frame = fits.getdata(self.data_dir + i)
            h = fits.open(self.data_dir +i)
            angle = h[0].header['PARANG'] + h[0].header['ROTPOSN'] - h[0].header['INSTANGL'] + zp_offset
            angle_list.append(angle)
        return angle_list

    #Takes an image and rotates it counterclockwise by the given parallactic angle 
    #https://github.com/AugustineStav/Angular-Differential-Imaging/blob/master/ADI_6-24-16.py

    def return_derotated_cube(self, cube, angle_list):

        '''
        This derotates the cube of images to get companion to match locations in every frame

        Args: 

        cube: cube of images (3d numpy array)
        angle_list: output from get_angle_list (list of angles to derotate each frame in the cube by)

        Returns:

        derotated_cube: a cube (3d numpy array) with de-rotation applied

        '''
        derotated_cube= []
        for i in range(cube.shape[0]):
            image = rotateImage(cube[i], -angle_list[i])
            derotated_cube.append(image)
            #print(centroid_1dg(image))
        return derotated_cube
            


    def shift_centroids(self, frames, target_location):

        '''
        This shifts the image using centroiding

        Args: 

        frames: 3d numpy array with frames to shift
        taret_location: where to place the star centroid. Tuple of (x pix, y pix)

        Returns:

        angle_list: python list with angles to distort each image in a cube

        '''
        centroids = []
        print(frames)
        for frame in frames:
            x,y = centroid_quadratic(frame)
            print(x,y)
            centroid = (x,y)
            centroids.append(centroid)
       # print(centroids)
        print(target_location)
        print(centroids)
        shifts = np.array(target_location) - np.array(centroids)
        shifted_frames = []
        for frame, shift_value in zip(frames, shifts):
            shifted_frame = scipy.ndimage.shift(frame, (shift_value[1],shift_value[0]))
            shifted_frames.append(shifted_frame)
        
        for i in shifted_frames:
            print(centroid_quadratic(i))

        return shifted_frames

    def run_it_all(self):

        '''
        This runs all the functions in sequence (i.e., runs the pipeline!)

        following the rxs example, after setting up the object simply run:

        rxs.run_it_all()
        
        '''


        sky, sky_median = self.get_cube_and_return_median(self.sky_filerange)
        target, target_median = self.get_cube_and_return_median(self.file_range)
        flats_on, flats_on_median = self.get_cube_and_return_median(self.flaton_filerange) #k band
        flats_off, flats_off_median = self.get_cube_and_return_median(self.flatoff_filerange)

        print("normalizing flat...")

        normalized_flat = (flats_on_median - flats_off_median)/np.median(flats_on_median - flats_off_median)

        print("finding bad pixel map...")

        index_flags = self.create_bad_pixel_map(normalized_flat)
        pixel_map = self.flag_pixels(index_flags, (1024, 1024))

        print("masking pixels from target, sky and flats...")

        target_masked = self.masked_cube(target, index_flags)
        flat_corrected = self.pixel_filter(normalized_flat, index_flags)
        sky_masked= self.pixel_filter(sky_median, index_flags)
        final_cube = (target_masked)/flat_corrected - sky_masked/flat_corrected


        print("writing final cube to fits files")

        os.chdir('/home/cdoo/nirc2_data/')
        pythoncode = os.system('mkdir ' + self.name)
        os.chdir('/home/cdoo/nirc2_data/' + self.name)

        for i in range(final_cube.shape[0]):
            hdu = fits.ImageHDU()
            hdu.data = np.array(final_cube[i])
            hdu.writeto(str(self.name) + "_" + str(i)+'.fits', overwrite = True)


        print('running distortion correction code...')

        rain_dir = self.distortion_correction()

        final_cube = self.get_final_cube(rain_dir)

        print('derotating images...')


        angle_list = self.get_angle_list()

        derotated_cube = self.return_derotated_cube(final_cube, angle_list)

        print('shifting centroids...')

        shifted_cube = self.shift_centroids(derotated_cube, (512, 512))

        median_cube = np.median(shifted_cube, axis = 0)

        hdu = fits.ImageHDU()
        hdu.data = np.array(median_cube)
        hdu.writeto(str(self.name) + '_final_median.fits', overwrite = 'True')


        hdu = fits.ImageHDU()
        hdu.data = np.array(shifted_cube)
        hdu.writeto(self.name + '_final_cube.fits', overwrite = 'True')


def rotateImage(image, angle):
    return skirotate(image, angle)


data_dir_feb = '/home/cdoo/nirc2_data/9feb/'
data_files_feb = np.array(sorted(os.listdir(data_dir_feb))) #get datafile list

rxs = Target("1rxs1216b", data_dir_feb, [79, 108], [109, 138], [3, 12], [13, 22])

rxs.run_it_all()


