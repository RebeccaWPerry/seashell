# start ipython with: ipython --gui=wx
# why the +20 and -20 offset between where the image is placed and where the points are placed?
# how to make things more transparent?

from mayavi import mlab
from numpy import ones, load, size
from numpy.random import rand

#a = batch_locate()
#numpy.save('shellstack_full_an1.npy',a)

#a = load('shellstack_full_an1')
b = load('shellstack_clustering.npy')

#mlab.figure(bgcolor=(1,1,1))
#numpoints = len(a) #34354
#mlab.points3d(a[:,0],a[:,1],a[:,2],4*ones(numpoints),resolution=16,scale_factor=1)

mlab.figure(bgcolor=(1,1,1))
numpoints = len(b) #16244
mlab.points3d(b[:,0],b[:,1],b[:,2],4*ones(numpoints),resolution=16,scale_factor=1,color=(.3,.3,0.3))



mlab.points3d(xs,ys,zs,4*ones(80),resolution=16,scale_factor=1,color=(1.0,0,.0))


#b = load('shellstack_clustering.npy')

#b = b[2000:11000]
#mlab.figure(bgcolor=(1,1,1))
#numpoints = len(b) #16244

'''final = load('outerclus.npy')
t = 0

for i in range(0,len(final)):
    cluster = final[i]
    numpoints = len(cluster)
    color = (rand(),rand(),rand())
    mlab.points3d(cluster[:,0],cluster[:,1],cluster[:,2],4*ones(numpoints),resolution=16,scale_factor=1,color=(1,0,0))
#mlab.points3d(b[:,0],b[:,1],b[:,2],4*ones(numpoints),resolution=16,scale_factor=1,color=(0,1,1))
'''
final = load('outerclus2.npy')
for i in range(0,len(final)):
    cluster = final[i]
    numpoints = len(cluster)
    color = (rand(),rand(),rand())
    mlab.points3d(cluster[:,0],cluster[:,1],cluster[:,2],4*ones(numpoints),resolution=16,scale_factor=1,color=color)
#mlab.points3d(b[:,0],b[:,1],b[:,2],4*ones(numpoints),resolution=16,scale_factor=1,color=(0,1,1))



#500 clusters


final = load('clusters2.npy')
t = 0

for i in range(0,len(final)):
    cluster = final[i]
    if size(cluster)>50:
        if (cluster[-1]-cluster[0])[2]>20:
	        print i
	        t+=1
	        numpoints = len(cluster)
	        color = (rand(),rand(),rand())
	        mlab.points3d(cluster[:,0],cluster[:,1],cluster[:,2],4*ones(numpoints),resolution=16,scale_factor=1,color=(1,0,0))


# one of the farthest out curves, plotted in black, 81 data points


ridge0 = final[789]
mlab.points3d(ridge0[:,0],ridge0[:,1],ridge0[:,2],4*ones(size(ridge0)/3),resolution=16,scale_factor=1,color=(0,0,0))


ridge1 = final[780]
mlab.points3d(ridge1[:,0],ridge1[:,1],ridge1[:,2],4*ones(size(ridge1)/3),resolution=16,scale_factor=1,color=(1,0,0))




filename = 'shellstack215.tif'
image1 = io.imread(filename)
image[image<12000]=False


obj = mlab.imshow(image1[0:,0:].T, opacity = .3, transparent = True)

#obj.actor.orientation = [0, 0, 0] # the required orientation 
obj.actor.position = [shape(image1)[0]/2.+20, shape(image1)[1]/2.-20, 215] # the required  position 
#obj.actor.scale = [0, 0, 0]


filename = 'shellstack345.tif'
image2 = io.imread(filename)
image[image<12000]=False


obj = mlab.imshow(image2[0:,0:].T, opacity = .3, transparent=True)

#obj.actor.orientation = [0, 0, 0] # the required orientation 
obj.actor.position = [shape(image2)[0]/2.+20, shape(image2)[1]/2.-20, 345] # the required  position 
#obj.actor.scale = [0, 0, 0]


stack = zeros([1375,1411,2])
stack[:,:,0]=image1
stack[:,:,1]=image2

#mlab.figure(bgcolor=(1, 1, 1))
#surf = mlab.surf(stack, colormap='cool')