 # Given a slice of a micro-CT scan of a shell, find the locations
 # of the bumps around the edge of the shell.
 # 
 # Given a whole stack of slices, locate the bumps in 3D by 
 # refining the bumps found in each slice.
 #
 # October, 2015
 # Rebecca W. Perry

 # TODO: make radius independent from padding
 # TODO: how do the results depend on choice of radius, searchpoints, 
        # and padding?
 # TODO: is oval the best shape? (see square example)
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

    [top, bottom] = [rowsum.min()-padding, rowsum.max()+padding+1]
    [left, right] = [colsum.min()-padding, colsum.max()+padding+1]
    
    assert top>=0, "padding value, %r pixels, is too large" % padding
    assert bottom<rows, "padding value, %r pixels, is too large" % padding
    assert left>=0, "padding value, %r pixels, is too large" % padding
    assert right<cols, "padding value, %r pixels, is too large" % padding

    smaller = binary_image[top:bottom, left:right]
    return smaller, [top, bottom, left, right]


def locate_bumps(image, radius=2, searchpoints=1000, display=False):
    '''Locates bumps around the exterior of a shape in a binarized image.
       Works in the frame of reference where the object is white and background 
       is black. Will convert images that don't follow this convention.

    Args:
        image (image or str): binary image in which to find bumps, or path to a 
        file containing the image (can be grayscale)

    Keyword arguments:
        radius (float): size of the oval around the object to calculate minimum 
        distances from. In units of 1/2 the length and width of cropped image.
        
        searchpoints (int): number of points to calculate distances from 
        (default 1000)
                
        display (bool): set to True to see plots of results (defalut False)

    Returns:
        2D ndarray listing coordinates (row, col, strength) of bumps'''
    #prepare image: load, threshold, crop, and find edges
    if type(image) == str:
        image = io.imread(image)
        binary = image > threshold_otsu(image)
    else:
        binary = image

    if binary[0,0] == True:
        binary = -binary

    binary, crop_location = crop(binary, padding = 2) #small pad for edges

    center = np.array(np.shape(binary))/2.
    edges = sobel(binary) #highlights pixels on either side of the edge
    edges = logical_and(edges,binary)
    edgepixs = np.array(np.where(edges>0)).T

    #find bumps using minimum distance to a bounding ellipse
    thetas = np.linspace(0,2*pi,searchpoints)
    ringpoints = np.array([center[0]+radius*center[0]*np.cos(thetas), 
                           center[1]+radius*center[1]*np.sin(thetas)]).T
    d = spatial.distance.cdist(ringpoints,edgepixs*1.0)
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
        final = zeros([len(set(part1)),3])

        for cluster in set(part1):
            final[cluster-1,0:2]  = mean(points[part1 == cluster],0)
            final[cluster-1,2] = sum(part1==cluster)

        ed = array(where(edges)).T
        t = spatial.distance.cdist(array([ed[:,0],ed[:,1]]).T,final[:,0:2])
        dist_from_border = amin(t,0)
        n = max(dist_from_border)
        num_clust1 += 1
    # end of reducing candidates


    # translate coordinates to their values in the whole frame
    master_coords = np.array([final[:,1]+crop_location[2],
                              final[:,0]+crop_location[0],final[:,2]]).T

    if display==True:

        plt.figure(figsize=[20,6])

        plt.subplot(1,4,1)
        io.imshow(image)
        plt.title('Raw Data')

        plt.subplot(1,4,2)
        io.imshow(binary)
        plt.plot(edgepixs[:,1],edgepixs[:,0],'g.')
        plt.plot(ringpoints[:,1],ringpoints[:,0],'g.')

        plt.subplot(1,4,3)
        io.imshow(binary)
        plt.gray()
        plt.plot(selected_candidates[1],selected_candidates[0],'ro')
        plt.xlim(0,np.shape(binary)[1])
        plt.ylim(np.shape(binary)[0],0)
        plt.title('First Guess at Peak Locations:\n'+str(np.shape(selected_candidates[1])[0])+ 
                  ' candidates' )

        plt.subplot(1,4,4)
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
        coords = locate_bumps(str(filenames[i]))
        if size(coords) > 0:
            this_framecoords3D = np.zeros([len(coords),4])
            this_framecoords3D[:,0:2] = coords[:,0:2]
            this_framecoords3D[:,2] = image_index
            this_framecoords3D[:,3] = coords[:,2]
            if size(coords3D)==0:
                coords3D = this_framecoords3D
            else:
                coords3D = np.append(coords3D, this_framecoords3D,0)
    return coords3D