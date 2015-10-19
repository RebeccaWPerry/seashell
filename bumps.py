 # Given a slice of a micro-CT scan of a shell, find the locations
 # of the bumps around the edge of the shell.
 # 
 # Given a whole stack of slices, locate the bumps in 3D by 
 # refining the bumps found in each slice.
 #
 # October, 2015
 # Rebecca W. Perry

 # TODO: make radius independent from padding
 # TODO: how do the results depend on choice of radius, searchpoints, and padding?
 # TODO: linking in 3D, more clustering algorithms?
 # TODO: read about convex hulls
 
from skimage import data, io
from skimage.filter import threshold_otsu, sobel
from scipy import spatial
from itertools import groupby
import matplotlib.pyplot as plt
from os import getcwd
import glob
from re import split
import numpy as np
import scipy.cluster.hierarchy as hac


def crop(binary_image, padding=0):
    """Crops the black or white border off of a binary image.

    Args: 
        binary_image (ndarray of bool): image to crop

    Keyword arguments:
        padding -- width of border to leave, in pixels (default 0)

    Returns:
        tuple including the subimage as an ndarray of size equal to or smaller 
        than the original binary_image and a list describing the range of the 
        subimage within the original binary image"""
    [rows, cols] = np.shape(binary_image)
    border_color = binary_image[0,0] #determines if the border is black or white
    rowsum = np.where(binary_image.sum(1)!=border_color*cols)[0]
    colsum = np.where(binary_image.sum(0)!=border_color*rows)[0]
    # TODO: return error if padding will cause indices outside image size
    [top, bottom] = [rowsum.min()-padding, rowsum.max()+padding+1]
    [left, right] = [colsum.min()-padding, colsum.max()+padding+1]
    smaller = binary_image[top:bottom, left:right]
    return smaller, [top, bottom, left, right]

def locate_bumps(image, radius=2, searchpoints=1000, display=False):
    '''Locates bumps around the exterior of a shape in an image.
    Args:
        image (str or binary image): image file or image in which to find bumps
        # TODO: allow for black on white or white on black
    Keyword arguments:
        radius (float): 
        
        searchpoints (int):
                
        display (bool): set to True to see plots of the results (defalut False)

    Returns:
        2D ndarray listing coordinates (row, col) of bumps'''
    #prepare image: load, threshold, crop, and find edges
    if type(image) == str:
        image = io.imread(image)
        binary = image > threshold_otsu(image)
    else:
        binary = image

    if sum(binary)/(1.0*np.size(binary))>.60:
        #sign that there is no object in the image
        master_coords= []
        print 'no object?'
    else:
        binary, crop_location = crop(binary, padding = 2)

        center = np.array(np.shape(binary))/2.
        edges = sobel(binary)
        edgepixs = np.array(np.where(edges>0)).T

        #find bumps using minimum distance to a bounding ellipse
        thetas = np.linspace(0,2*pi,searchpoints)
        ringpoints = np.array([center[1]+radius*center[1]*np.cos(thetas), 
                               center[0]+radius*center[0]*np.sin(thetas)]).T
        d = spatial.distance.cdist(ringpoints,edgepixs)
        votes = edgepixs[np.argmin(d,1)] #each ringpoint votes once


        # group candidates that are very close together, and probably represent the same bump
        selected_candidates = votes.T
        points = votes
        z = hac.linkage(points, method='single')
        knee = -1 * np.diff(z[::-1, 2], 1)
        num_clust1 = knee[10:].argmax() + 10 + 1 # guess of number of clusters

        n = 3
        while n>2: #don't let any center of mass be more than 2 pixels from the original edge
            part1 = hac.fcluster(z, num_clust1, 'maxclust')
            final = zeros([len(set(part1)),2])

            for cluster in set(part1):
                final[cluster-1]  = mean(points[part1 == cluster],0)

            ed = array(where(edges)).T
            t = spatial.distance.cdist(array([ed[:,0],ed[:,1]]).T,final)
            dist_from_border = amin(t,0)
            print dist_from_border
            n = max(dist_from_border)
            num_clust1 += 1
        # end of reducing candidates


        # reset coordinates to the whole frame
        master_coords = np.array([final[:,1]+crop_location[2],final[:,0]+crop_location[0]]).T

        if display==True:

            plt.figure(figsize=[20,6])

            plt.subplot(1,3,1)
            io.imshow(image)
            plt.title('Raw Data')

            plt.subplot(1,3,2)
            io.imshow(binary)
            plt.gray()
            plt.plot(selected_candidates[1],selected_candidates[0],'ro')
            plt.xlim(0,np.shape(binary)[1])
            plt.ylim(np.shape(binary)[0],0)
            plt.title('First Guess at Peak Locations:\n'+str(np.shape(selected_candidates[1])[0])+ 
                      ' candidates' )

            plt.subplot(1,3,3)
            io.imshow(image)
            plt.plot(master_coords[:,0],master_coords[:,1],'ro')
            plt.title('Refined Peak Locations Shown on Raw Data:\n'+str(np.shape(final)[0])+ 
                      ' peaks')
            plt.xlim(0,np.shape(image)[1])
            plt.ylim(np.shape(image)[0],0)

            io.show()

    return master_coords

def batch_locate(directory='here'):
    '''Batch processin version of locate_bumps.
    Keyword arguments:
        directory (str): directory in which .tif images are stored. Default
            is current working directory.

    Returns:
        2D ndarray listing coordinates (row, col, image number) of bumps'''
    # glob file names, sort
    if directory == 'here':
        directory = getcwd()
    filenames = np.sort(glob.glob(directory+'/*.tif'))

    if len(filenames) == 0:
        raise SystemExit('No tif files were found in the specified directory.')

    coords3D = []

    for i in np.arange(0,len(filenames)):
        image_index = int(split('(\d+)',filenames[i].split(directory)[1])[1])
        print ('Image Number: ' + str(image_index))
        coords = locate_bumps(filenames[i])
        if size(coords) > 0:
            this_framecoords3D = np.zeros([len(coords),3])
            this_framecoords3D[:,0:2] = coords
            this_framecoords3D[:,2] = image_index
            if size(coords3D)==0:
                coords3D = this_framecoords3D
            else:
                coords3D = np.append(coords3D, this_framecoords3D,0)
    return coords3D